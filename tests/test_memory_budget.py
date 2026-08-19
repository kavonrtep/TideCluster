#!/usr/bin/env python3
"""Unit tests for the memory budget resolution chain (issue #6).

The pool-sizing gates (TideHunter worker pool, TAREAN thread cap) used to read
``/proc/meminfo`` ``MemAvailable`` as their only budget source. ``MemAvailable``
is not namespaced by the kernel, so under a batch scheduler or inside a
container it reports the *host's* free memory: a 128 GB PBS job saw a 1.6 TB
budget, ran 32 x ~9.7 GB TideHunter parts concurrently and was OOM-killed hours
in, with nothing warning beforehand.

These tests pin the parts that are easy to get silently wrong:
  1. precedence: --max_memory > AGENT_MEMORY > scheduler env > cgroup > MemAvailable
  2. unit conversion per source (bytes vs KB vs MB vs GB) -- a 1024x slip here
     re-creates the original bug
  3. "no limit" sentinels: cgroup v2 "max" and the v1 huge-number sentinel must
     both read as *absent*, not as an enormous budget
  4. the cgroup chain: a limit on a parent slice must be found, and the
     effective (smallest) limit wins
  5. the scheduler + MemAvailable warning that makes this diagnosable

Run: python3 tests/test_memory_budget.py
"""
import os
import shutil
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import tc_utils as tc

H = tc.MEMORY_HEADROOM
GB = 1024.0
NOWHERE = "/nonexistent-tidecluster-test-path"

failures = []


def check(label, got, want, tol=1e-6):
    ok = (got == want if not isinstance(want, float)
          else abs(got - want) <= tol * max(1.0, abs(want)))
    print(F"  {'ok  ' if ok else 'FAIL'} {label}: got {got!r}, want {want!r}")
    if not ok:
        failures.append(label)


def budget(env, sysfs=NOWHERE, meminfo=NOWHERE, proc_cgroup=NOWHERE,
           explicit=None):
    """memory_budget_mb with every filesystem source disabled unless given."""
    return tc.memory_budget_mb(explicit_gb=explicit, environ=env,
                               sysfs_root=sysfs, proc_cgroup=proc_cgroup,
                               meminfo=meminfo)


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def fake_meminfo(tmp, kb):
    path = os.path.join(tmp, "meminfo")
    write(path, F"MemTotal:       {kb * 2} kB\nMemAvailable:   {kb} kB\n"
                "SwapTotal:             0 kB\n")
    return path


# --- 1. precedence -----------------------------------------------------------
print("precedence (first hit wins)")
tmp = tempfile.mkdtemp(prefix="tc_membudget_")
try:
    meminfo = fake_meminfo(tmp, 200 * 1024 * 1024)          # 200 GB available
    sysfs = os.path.join(tmp, "cgroup_v2")
    write(os.path.join(sysfs, "memory.max"), str(64 * 1024 * 1024 * 1024))

    # every source present at once, peeled off one at a time
    full = {"AGENT_MEMORY": "48", "PBS_RESC_MEM": str(32 * 1024 ** 3),
            "SLURM_MEM_PER_NODE": "16384", "LSB_MAX_MEM_RUSAGE": "8388608"}

    mb, src = budget(full, sysfs=sysfs, meminfo=meminfo, explicit=128)
    check("--max_memory beats everything", (mb, src), (128 * GB * H, "--max_memory"))

    mb, src = budget(full, sysfs=sysfs, meminfo=meminfo)
    check("AGENT_MEMORY beats scheduler env", (mb, src), (48 * GB * H, "AGENT_MEMORY"))

    env = dict(full); env.pop("AGENT_MEMORY")
    mb, src = budget(env, sysfs=sysfs, meminfo=meminfo)
    check("PBS_RESC_MEM beats other scheduler vars",
          (mb, src), (32 * GB * H, "PBS_RESC_MEM"))

    env.pop("PBS_RESC_MEM")
    mb, src = budget(env, sysfs=sysfs, meminfo=meminfo)
    check("SLURM_MEM_PER_NODE beats LSF", (mb, src), (16 * GB * H, "SLURM_MEM_PER_NODE"))

    env.pop("SLURM_MEM_PER_NODE")
    mb, src = budget(env, sysfs=sysfs, meminfo=meminfo)
    check("LSB_MAX_MEM_RUSAGE beats cgroup", (mb, src), (8 * GB * H, "LSB_MAX_MEM_RUSAGE"))

    mb, src = budget({}, sysfs=sysfs, meminfo=meminfo)
    check("cgroup beats MemAvailable", (mb, src), (64 * GB * H, "cgroup"))

    mb, src = budget({}, meminfo=meminfo)
    check("MemAvailable is the last resort", (mb, src), (200 * GB * H, "MemAvailable"))

    mb, src = budget({})
    check("nothing readable -> (None, 'none')", (mb, src), (None, "none"))

    # --- 2. unit conversions -------------------------------------------------
    print("unit conversions (a 1024x slip re-creates the bug)")
    check("PBS_RESC_MEM is bytes",
          budget({"PBS_RESC_MEM": str(128 * 1024 ** 3)})[0], 128 * GB * H)
    check("SLURM_MEM_PER_NODE is MB",
          budget({"SLURM_MEM_PER_NODE": "131072"})[0], 128 * GB * H)
    check("LSB_MAX_MEM_RUSAGE is KB",
          budget({"LSB_MAX_MEM_RUSAGE": str(128 * 1024 ** 2)})[0], 128 * GB * H)
    check("AGENT_MEMORY is GB", budget({"AGENT_MEMORY": "128"})[0], 128 * GB * H)
    check("SLURM_MEM_PER_CPU x SLURM_CPUS_ON_NODE",
          budget({"SLURM_MEM_PER_CPU": "4096", "SLURM_CPUS_ON_NODE": "32"}),
          (128 * GB * H, "SLURM_MEM_PER_CPU"))
    check("explicit unit suffix honoured over the default unit",
          budget({"PBS_RESC_MEM": "128gb"})[0], 128 * GB * H)
    check("MemAvailable is KB",
          budget({}, meminfo=fake_meminfo(tmp, 128 * 1024 * 1024))[0], 128 * GB * H)

    print("malformed values fall through instead of aborting")
    check("non-numeric AGENT_MEMORY skipped",
          budget({"AGENT_MEMORY": "lots", "SLURM_MEM_PER_NODE": "16384"}),
          (16 * GB * H, "SLURM_MEM_PER_NODE"))
    check("empty scheduler var skipped", budget({"PBS_RESC_MEM": "  "}), (None, "none"))
    check("zero scheduler var skipped", budget({"PBS_RESC_MEM": "0"}), (None, "none"))
    check("SLURM_MEM_PER_CPU without CPUS_ON_NODE skipped",
          budget({"SLURM_MEM_PER_CPU": "4096"}), (None, "none"))

    # --- 3. cgroup sentinels ------------------------------------------------
    print("cgroup 'no limit' sentinels must read as absent")
    v2max = os.path.join(tmp, "cg_v2_max")
    write(os.path.join(v2max, "memory.max"), "max\n")
    check("v2 'max' -> no limit",
          tc._cgroup_limit_mb(sysfs_root=v2max, proc_cgroup=NOWHERE), None)
    check("v2 'max' falls through to MemAvailable",
          budget({}, sysfs=v2max, meminfo=fake_meminfo(tmp, 10 * 1024 * 1024))[1],
          "MemAvailable")

    v1inf = os.path.join(tmp, "cg_v1_inf")
    write(os.path.join(v1inf, "memory", "memory.limit_in_bytes"),
          "9223372036854771712\n")            # v1 "unlimited" sentinel
    check("v1 huge sentinel -> no limit",
          tc._cgroup_limit_mb(sysfs_root=v1inf, proc_cgroup=NOWHERE), None)

    v1 = os.path.join(tmp, "cg_v1")
    write(os.path.join(v1, "memory", "memory.limit_in_bytes"),
          str(96 * 1024 ** 3))
    check("v1 finite limit read",
          tc._cgroup_limit_mb(sysfs_root=v1, proc_cgroup=NOWHERE), 96 * GB)

    # --- 4. cgroup chain ----------------------------------------------------
    print("cgroup chain: parent-slice limit found, smallest wins")
    chain = os.path.join(tmp, "cg_chain")
    proc = os.path.join(tmp, "proc_cgroup")
    write(proc, "0::/user.slice/job-42.scope\n")
    write(os.path.join(chain, "memory.max"), "max\n")                     # root
    write(os.path.join(chain, "user.slice", "memory.max"),
          str(128 * 1024 ** 3))                                           # parent
    os.makedirs(os.path.join(chain, "user.slice", "job-42.scope"))
    write(os.path.join(chain, "user.slice", "job-42.scope", "memory.max"),
          "max\n")                                                        # leaf
    check("limit on a parent slice is found",
          tc._cgroup_limit_mb(sysfs_root=chain, proc_cgroup=proc), 128 * GB)

    write(os.path.join(chain, "user.slice", "job-42.scope", "memory.max"),
          str(64 * 1024 ** 3))
    check("smallest limit in the chain wins",
          tc._cgroup_limit_mb(sysfs_root=chain, proc_cgroup=proc), 64 * GB)

    # A container's namespaced root: /proc/self/cgroup names a path that does
    # not resolve under /sys/fs/cgroup, but the root itself carries the limit.
    ns = os.path.join(tmp, "cg_ns")
    ns_proc = os.path.join(tmp, "proc_cgroup_ns")
    write(ns_proc, "0::/user.slice/user-1000.slice/apptainer-2026625.scope\n")
    write(os.path.join(ns, "memory.max"), str(128 * 1024 ** 3))
    check("namespaced root limit found when the path does not resolve",
          tc._cgroup_limit_mb(sysfs_root=ns, proc_cgroup=ns_proc), 128 * GB)

    print("cgroup v1 multi-controller /proc/self/cgroup line")
    v1chain = os.path.join(tmp, "cg_v1_chain")
    v1proc = os.path.join(tmp, "proc_cgroup_v1")
    write(v1proc, "8:cpu,cpuacct:/torque/999.pbs\n"
                  "7:memory:/torque/999.pbs\n"
                  "0::/\n")
    write(os.path.join(v1chain, "memory", "memory.limit_in_bytes"),
          "9223372036854771712\n")                                        # root
    write(os.path.join(v1chain, "memory", "torque", "999.pbs",
                       "memory.limit_in_bytes"), str(128 * 1024 ** 3))
    check("v1 leaf limit found via the memory controller path",
          tc._cgroup_limit_mb(sysfs_root=v1chain, proc_cgroup=v1proc), 128 * GB)

    # --- 5. the warning -----------------------------------------------------
    print("host-budget warning under a scheduler")
    tc._warned_host_budget = False
    check("MemAvailable + PBS_JOBID warns",
          tc.warn_if_host_memory_budget("MemAvailable", {"PBS_JOBID": "42.pbs"}), True)
    check("warns only once per process",
          tc.warn_if_host_memory_budget("MemAvailable", {"PBS_JOBID": "42.pbs"}), False)
    tc._warned_host_budget = False
    check("MemAvailable with no scheduler stays quiet",
          tc.warn_if_host_memory_budget("MemAvailable", {}), False)
    check("a real budget source stays quiet",
          tc.warn_if_host_memory_budget("cgroup", {"SLURM_JOB_ID": "7"}), False)
    check("--max_memory stays quiet",
          tc.warn_if_host_memory_budget("--max_memory", {"PBS_JOBID": "42.pbs"}), False)
    tc._warned_host_budget = False

    # --- 6. what the gates do with it ---------------------------------------
    # The whole point: the PBS case from issue #6 must pick a pool that fits.
    print("issue #6 regression: 128 GB job, 9.7 GB/part, 32 cores")
    mb, src = budget({"PBS_RESC_MEM": str(128 * 1024 ** 3), "PBS_JOBID": "1.pbs"},
                     meminfo=fake_meminfo(tmp, 1600 * 1024 * 1024))
    pool = max(1, min(32, int(mb // 9736)))
    check("budget source is the job's, not the host's", src, "PBS_RESC_MEM")
    check("pool_size fits under the cgroup limit", pool, 10)
    check("pool_size x per-part peak <= real limit", pool * 9736 <= 128 * GB, True)

    # --- 7. the TAREAN thread cap uses the same chain ------------------------
    # 32 x 4 GB = 128 GB would have been the next thing to OOM in issue #6.
    print("TAREAN thread cap (TideCluster._tarean_max_threads)")
    sys.argv = [sys.argv[0]]          # TideCluster.py parses argv only under main
    import TideCluster as tcl
    saved_env = os.environ.copy()
    try:
        os.environ.clear()
        os.environ.update({"PBS_RESC_MEM": str(128 * 1024 ** 3),
                           "PBS_JOBID": "1.pbs"})
        threads, mb, src = tcl._tarean_max_threads(32)
        check("cap comes from the job's limit", src, "PBS_RESC_MEM")
        check("threads capped below the core count", threads, 26)
        check("threads x 4 GB <= real limit", threads * 4000 <= 128 * GB, True)

        os.environ.clear()
        threads, mb, src = tcl._tarean_max_threads(32, max_memory=16)
        check("--max_memory caps threads", (threads, src), (3, "--max_memory"))

        os.environ.clear()
        os.environ["AGENT_MEMORY"] = "1024"
        threads, mb, src = tcl._tarean_max_threads(8)
        check("a large budget leaves -c as the only cap", threads, 8)
    finally:
        os.environ.clear()
        os.environ.update(saved_env)
finally:
    shutil.rmtree(tmp, ignore_errors=True)

print()
if failures:
    print(F"FAILED ({len(failures)}): " + ", ".join(failures))
    sys.exit(1)
print("test_memory_budget.py: all checks passed")

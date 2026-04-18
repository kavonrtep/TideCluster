#!/bin/sh
# Install tidecluster into $PREFIX. noarch package; copies the tree into
# $PREFIX/share/tidecluster and symlinks every CLI into $PREFIX/bin.
# Mirrors the symlink set in kavonrtep/recipes/tidecluster/build.sh.
set -x -e

PKG_DIR="${PREFIX}/share/tidecluster"
mkdir -p "${PREFIX}/bin" "${PKG_DIR}"
cp -r . "${PKG_DIR}"

ln -s "${PKG_DIR}/TideCluster.py"                       "${PREFIX}/bin/TideCluster.py"
ln -s "${PKG_DIR}/tc_update_gff3.py"                    "${PREFIX}/bin/tc_update_gff3.py"
ln -s "${PKG_DIR}/tc_reannotate.py"                     "${PREFIX}/bin/tc_reannotate.py"
ln -s "${PKG_DIR}/tc_merge_annotations.py"              "${PREFIX}/bin/tc_merge_annotations.py"
ln -s "${PKG_DIR}/tc_utils.py"                          "${PREFIX}/bin/tc_utils.py"
ln -s "${PKG_DIR}/tc_comparative_analysis.R"            "${PREFIX}/bin/tc_comparative_analysis.R"
ln -s "${PKG_DIR}/tc_summarize_comparative_analysis.R"  "${PREFIX}/bin/tc_summarize_comparative_analysis.R"

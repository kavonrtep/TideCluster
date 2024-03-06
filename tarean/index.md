<style>
     p {
         max-width: 700px;
         text-align: left;
     }
     ul {
         max-width: 700px;
         text-align: left;
     }
     ul li {
         margin-bottom: 0px;
         padding: 5px;
     }   

     body {
       font-family: Verdana, sans-serif;
       font-size: 13px;
       width: 100%;
       max-width: 700px;
       margin: 30px;
       display: block;
       text-align: left;
     }
     .frame {

         width: 100%;
         border: 0px;
         white-space: nowrap;
         text-align: left; 
         margin: 1em 0;
     }
     .helper {
         display: inline-block;
         height: 100%;
         vertical-align: top;
     }
     img {
         background: #FFFFFF;
         vertical-align: top;
     }

     table {
         border-collapse: collapse;
         width: 100%;
     }
     th.date {
         width: 20%
     }
     th.subject {
         width: 70%; /* Not necessary, since only 70% width remains */
     }

     tr {
         border-bottom: 1px solid #ccc;
     }

     th {
         text-align: left;    
     }
     th, td {
         padding: 5px;
     }
</style>
# TideCluster report

## Analysis settings:

    SETTINGS_PLACEHOLDER

## Summary

    SUMMARY_PLACEHOLDER

## Tandem Repeat Analysis (TAREAN)

[**Report**](PREFIX_PLACEHOLDER_tarean_report.html)

Report includes TAREAN analysis for all tandem repeat clusters (TRCs) pass the threshold of the minimamal combined length of tandem repeat arrays within a single cluster

- **Consensus sequence** for each tandem repeat cluster (TRC) estimated by k-mer-based de Brujn graph approach.
- **Consensus sequence for simple sequence repeats** (SSRs) estimated by simplified k-mer-based approach.
- Basic **TRC statistics** (number of arrays, total length, minimal and maximal length of arrays)


## K-Mer Interval Tandem Repeat Estimation (KITE)

[**Report**](PREFIX_PLACEHOLDER_kite_report.html)

The K-Mer Interval Tandem Repeat Estimation (KITE) method is applied to analyze individual Tandem Repeat Arrays (TRAs) within each TRC. In contrast to TAREAN, which provides monomer size estimates for the entire TRC, KITE also estimates monomer sizes for each array by evaluating the distances between all repetitive k-mers across the tandem repeat array, with each tandem array analyzed individually. This method facilitates the detection of higher-order repeats and captures the variability in monomer size across different tandem repeat arrays. The KITE analysis includes: 

- **Monomer size estimate for the entire TRC** : This estimate may differ from the TAREAN estimate due to the different methodologies employed. 
- **Monomer size estimates for individual TRAs** : These estimates are provided as: 
  - **Primary estimate** : The estimate with the highest score. 
  - **Alternative estimates** : Estimates with lower scores, these estimates help in identifying **higher-order repeats**.



## Credits

TideCluster utilizes Tidehunter [https://github.com/Xinglab/TideHunter] for tandem repeat detection and TAREAN for reconstruction of consensus sequences of tandem repeats.
If you use TideCluster please cite:

- https://github.com/kavonrtep/TideCluster [![DOI](https://zenodo.org/badge/601111441.svg)](https://zenodo.org/badge/latestdoi/601111441)
- TAREAN: a computational tool for identification and characterization of satellite DNA from unassembled short reads (https://doi.org/10.1093/nar/gkx257) 
- TideHunter: efficient and sensitive tandem repeat detection from noisy long-reads using seed-and-chain (https://doi.org/10.1093/bioinformatics/btz376)

### Authors

Petr Novak, Jiri Macas,  Laboratory of Molecular Cytogenetics, Biology Centre CAS, Czech Republic

# link of websummary metrics intepration
https://www.10xgenomics.com/support/universal-three-prime-gene-expression/documentation/steps/sequencing/interpreting-cell-ranger-web-summary-files-for-single-cell-gene-expression-assays

1. Important metrics on measure overall success
Several Metrics in the web summary file can be used to assess the overall success of an experiment, 
including sequencing, mapping, and cell metrics.

# First pay attention on the color of Main numers and the "Alerts" section on the top
Color: green text indicates that key metrics are in the expected range
red/yellow: indicateds errors/warnings
Descriptions of the metrics can also be found by clicking the icon ? next to the section header.

# sequencing metrics: measure the sequencing peformance
total reads: depending on the sequencing depth
valid barcodes: > 75%
valid UMI: >75%
sequencing saturation: dpendent  on sequence depth and sample complexity
Q30 bases in RNA reads: sequencing platform dependent, but ideally >65%

# cell metrics: measure the sample quality
** estimated cell numbers: depending on the input cells in the exp design, 500-10,000, 
  higher or lower than expected indicating inaccurate cell counting, cell lysis, or failure during GEM generation
mean reads per cell: 

** fraction reads in cells: >70%
  The fraction of reads that contain a valid barcode and confidenely mapped to the transcriptiome and the associated barcode is called as a cell
  lower percentages indicate a high level of ambient RNA partitioned into all GEMs

** median reads per cell: recommend minimial 20K reads per cell
  Here our samples are 45K reads per cell
  
** median genes per cell: depending on the cell type and sequencing depth
  ~2000 is OK

** total genes detected: depending on the cell type and sequencing depth
  20k-24k is OK

# mapping metrics
Reads mapped confidently to transcriptome: variable but ideally > 30%
Reads mapped antisense to gene: ideal < 10%, the value may be higher if using a pre-mRNA reference or may indicate incorrect Gel Bead chemistry
Other metrics: variable

2. Check the plots
# barcode rank plot:
looking for a steep drop-off (an indication of good seperation between cell-associated barcodes and barcodes within empty GEMs)
if biomodal (having two cliff and knee) distribution, indiates the heterogenous population of cells in a sample
if round curve and llack of steep cliff, indicate low sample quality or loss of single-cell behavior, e.g., wetting failure, premature cell lysis, low cell viability
if having a steep drop-off but the total number of barcode signiffcantly lower than expected, indicating a sample clog or inaccurate cell counting


# t-SNE plots:
look for strucutred clusters with clear sepration between high UMI and low UMI containing barcodes


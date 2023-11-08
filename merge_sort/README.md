# Merge Sort

This approach divides the array of numbers based on the number of processors available. Then, merge sort is carried out on these sub-arrays. When these are sorted, the original array is merged back together.

The following times were observed per rank, when sorting 10,000,000 elements. Execution times is split between communication operations and computation operations.

### 1 Processor
|Processor|Comm. Time (s)|Comp. Time (s)|Total Time (s)|
|---|----|---|---|
|0|0.042|2.772|2.814|

### 2 Processors
|Processor|Comm. Time (s)|Comp. Time (s)|Total Time (s)|
|---|----|---|---|
|Total|0.264|2.797|3.061|
|Average|0.132|1.4|1.532|
|0|0.037|1.463|1.5|
|1|0.227|1.334|1.561|

### 4 Processors
|Processor|Comm. Time (s)|Comp. Time (s)|Total Time (s)|
|---|----|---|---|
|Total|0.867|3.052|3.919|
|Average|0.217|0.763|0.980|
|0|0.039|0.994|1.038|
|1|0.287|0.668|0.955|
|2|0.269|0.690|0.959|
|3|0.272|0.700|1.266|
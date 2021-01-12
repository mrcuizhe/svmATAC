# svmATAC
Enhancement and integration of peak signal enables accurate identification of cell type in scATAC-seq

## Dependencies

- R 3.6
- Python 3.7

## Quick Start

The pipeline of svmATAC is consist of two steps:

- 1.[Pre-process](#Pre-process):
	- If you are trying to reproduce the results in svmATAC paper from **raw sequencing data**(Can be found in preprocess/XXX/input folder), you may need this step. 
	- If you are trying to reproduce the results in svmATAC paper using **our provided data**(Can be found in XXX-dataset/XXX/input folder), or your data are already merged to 0-1 matrix (peak-cell) and labelled, you may skip this step.

- 2.[Training-Classification](#Training-Classification):
	- The training scripts and classification scripts are stored in bin/ folder of each experiment, you can try these scripts and reproduce the results in svmATAC paper through following chapters: 
		- [intra-dataset experiment](#intra-dataset)
			- Corces2016
			- Buenrostro2018
			- 10xPBMCsV1
			- 10xPBMCsNextGem
			- 10xPBMCsV1-labeled
			- 10xPBMCsNextGem-labeled
		- [inter-dataset experiment](#inter-dataset)
			- 10xPBMCs-labeled
			- 10xPBMCs-unlabeled
	- For each single experiment:
		- All the scripts are available in the bin/ folder
			- All scripts are numbered and users should execute one by one.
			- The intermediate temporary files are stored in tmp/ folder and you can ignore these files.
		- All the input data required are available in the input/ folder
		- All the output data generated are stored in the output/ folder

## Content	
	
### Pre-process

This chapter stores the scripts for processing raw data.

- The _**Corces2016**_ dataset (folder **'Corces2016'**)
	- [bin](https://github.com/mrcuizhe/svmATAC/tree/master/preprocess/Corces2016/bin)
	- [input](https://github.com/mrcuizhe/svmATAC/tree/master/preprocess/Corces2016/input)
	- [output](https://github.com/mrcuizhe/svmATAC/tree/master/preprocess/Corces2016/output)

- The _**Buenrostro2018**_ dataset (folder **'Buenrostro2018'**)
	- [bin](https://github.com/mrcuizhe/svmATAC/tree/master/preprocess/Buenrostro2018/bin)
	- [input](https://github.com/mrcuizhe/svmATAC/tree/master/preprocess/Buenrostro2018/input)
	- [output](https://github.com/mrcuizhe/svmATAC/tree/master/preprocess/Buenrostro2018/output)

- The _**10xPBMCsV1**_ dataset (folder **'10xPBMCsV1'**)
	- [bin](https://github.com/mrcuizhe/svmATAC/tree/master/preprocess/10xPBMCsV1/bin)
	- [input](https://github.com/mrcuizhe/svmATAC/tree/master/preprocess/10xPBMCsV1/input)
	- [output](https://github.com/mrcuizhe/svmATAC/tree/master/preprocess/10xPBMCsV1/output)

- The _**10xPBMCsNextGem**_ dataset (folder **'10xPBMCsNextGem'**)
	- [bin](https://github.com/mrcuizhe/svmATAC/tree/master/preprocess/10xPBMCsNextGem/bin)
	- [input](https://github.com/mrcuizhe/svmATAC/tree/master/preprocess/10xPBMCsNextGem/input)
	- [output](https://github.com/mrcuizhe/svmATAC/tree/master/preprocess/10xPBMCsNextGem/output)

This chapter stores the scripts for assigning labels to 10x PBMCs v1 and nextGem scATAC-seq data from labeled scRNA-seq data using Seurat.

- The _**scRNA-seq-5k-v3**_ dataset (folder **'scRNA-seq-5k-v3'**)
	- [bin](https://github.com/mrcuizhe/svmATAC/tree/master/10x_assign_label/scRNA-seq-5k-v3/bin)
	- [data](https://github.com/mrcuizhe/svmATAC/tree/master/10x_assign_label/scRNA-seq-5k-v3/data)
	- [output](https://github.com/mrcuizhe/svmATAC/tree/master/10x_assign_label/scRNA-seq-5k-v3/output)

- The _**scATAC-seq-5k-v1**_ dataset (folder **'scATAC-seq-5k-v1'**)
	- [bin](https://github.com/mrcuizhe/svmATAC/tree/master/10x_assign_label/scATAC-seq-5k-v1/bin)
	- [input](https://github.com/mrcuizhe/svmATAC/tree/master/10x_assign_label/scATAC-seq-5k-v1/input)
	- [output](https://github.com/mrcuizhe/svmATAC/tree/master/10x_assign_label/scATAC-seq-5k-v1/output)

- The _**scATAC-seq-5k-nextgem**_ dataset (folder **'scATAC-seq-5k-nextgem'**)
	- [bin](https://github.com/mrcuizhe/svmATAC/tree/master/10x_assign_label/scATAC-seq-5k-nextgem/bin)
	- [input](https://github.com/mrcuizhe/svmATAC/tree/master/10x_assign_label/scATAC-seq-5k-nextgem/input)
	- [output](https://github.com/mrcuizhe/svmATAC/tree/master/10x_assign_label/scATAC-seq-5k-nextgem/output)

### Training-Classification

#### intra-dataset

This chapter stores the scripts for intra-dataset experiments which are described in manuscript.

- The _**Corces2016**_ dataset (folder **'Corces2016'**)
	- [bin](https://github.com/mrcuizhe/svmATAC/tree/master/intra-dataset/Corces2016/bin)
	- [input](https://github.com/mrcuizhe/svmATAC/tree/master/intra-dataset/Corces2016/input)
	- [output](https://github.com/mrcuizhe/svmATAC/tree/master/intra-dataset/Corces2016/output)

- The _**Buenrostro2018**_ dataset (folder **'Buenrostro2018'**)
	- [bin](https://github.com/mrcuizhe/svmATAC/tree/master/intra-dataset/Buenrostro2018/bin)
	- [input](https://github.com/mrcuizhe/svmATAC/tree/master/intra-dataset/Buenrostro2018/input)
	- [output](https://github.com/mrcuizhe/svmATAC/tree/master/intra-dataset/Buenrostro2018/output)

- The _**10xPBMCsV1**_ dataset (folder **'10xPBMCsV1'**)
	- [bin](https://github.com/mrcuizhe/svmATAC/tree/master/intra-dataset/10xPBMCsV1/bin)
	- [input](https://github.com/mrcuizhe/svmATAC/tree/master/intra-dataset/10xPBMCsV1/input)
	- [output](https://github.com/mrcuizhe/svmATAC/tree/master/intra-dataset/10xPBMCsV1/output)

- The _**10xPBMCsNextGem**_ dataset (folder **'10xPBMCsNextGem'**)
	- [bin](https://github.com/mrcuizhe/svmATAC/tree/master/intra-dataset/10xPBMCsNextGem/bin)
	- [input](https://github.com/mrcuizhe/svmATAC/tree/master/intra-dataset/10xPBMCsNextGem/input)
	- [output](https://github.com/mrcuizhe/svmATAC/tree/master/intra-dataset/10xPBMCsNextGem/output)
	
- The _**10xPBMCsV1-labeled**_ dataset (folder **'10xPBMCsV1_labeled'**)
	- [bin](https://github.com/mrcuizhe/svmATAC/tree/master/intra-dataset/10xPBMCsV1_labeled/bin)
	- [input](https://github.com/mrcuizhe/svmATAC/tree/master/intra-dataset/10xPBMCsV1_labeled/input)
	- [output](https://github.com/mrcuizhe/svmATAC/tree/master/intra-dataset/10xPBMCsV1_labeled/output)

- The _**10xPBMCsNextGem-labeled**_ dataset (folder **'10xPBMCsNextGem_labeled'**)
	- [bin](https://github.com/mrcuizhe/svmATAC/tree/master/intra-dataset/10xPBMCsNextGem_labeled/bin)
	- [input](https://github.com/mrcuizhe/svmATAC/tree/master/intra-dataset/10xPBMCsNextGem_labeled/input)
	- [output](https://github.com/mrcuizhe/svmATAC/tree/master/intra-dataset/10xPBMCsNextGem_labeled/output)
	
#### inter-dataset

This chapter stores the scripts for inter-dataset experiments which are described in manuscript.

- The _**10xPBMCs-labeled**_ dataset (folder **'10xPBMCs_labeled'**)
	- [bin](https://github.com/mrcuizhe/svmATAC/tree/master/inter-dataset/10xPBMCs_labeled/bin)
	- [input](https://github.com/mrcuizhe/svmATAC/tree/master/inter-dataset/10xPBMCs_labeled/input)
	- [output](https://github.com/mrcuizhe/svmATAC/tree/master/inter-dataset/10xPBMCs_labeled/output)
	
- The _**10xPBMCs-unlabeled**_ dataset (folder **'10xPBMCs_unlabeled'**)
	- [bin](https://github.com/mrcuizhe/svmATAC/tree/master/inter-dataset/10xPBMCs_unlabeled/bin)
	- [input](https://github.com/mrcuizhe/svmATAC/tree/master/inter-dataset/10xPBMCs_unlabeled/input)
	- [output](https://github.com/mrcuizhe/svmATAC/tree/master/inter-dataset/10xPBMCs_unlabeled/output)
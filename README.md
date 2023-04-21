# AWS Batch Architecture for Protein Folding and Design

## Contents

1. [Overview](#1-overview)  
2. [Quick Start](#2-quick-start)  
3. [Advanced Configuration](#3-advanced-configuration)  
    3.1. [Optional CloudFormation Parameters](#31-optional-cloudformation-parameters)  
    3.2. [Manual Data Download](#32-manual-data-download)  
    3.3. [Clean Up](#33-clean-up)  
4. [Module Information](#4-module-information)  
    4.1. [JackHMMER](#41-jackhmmer)  
    4.2. [AlphaFold](#42-alphafold)  
    4.3. [OpenFold](#43-openfold)  
    4.4. [OmegaFold](#44-omegafold)  
    4.5. [RFDesign](#45-rfdesign)  
    4.6. [ESMFold](#46-esmfold)  
    4.7. [ProteinMPNN](#47-proteinmpnn)  
    4.8. [DiffDock](#48-diffdock)  
    4.9. [RFDiffusion](#49-rfdiffusion)  
5. [Architecture Details](#5-architecture-details)  
    5.1. [Stack Creation Details](#51-stack-creation-details)  
    5.2. [Cost](#52-cost)  
6. [FAQ](#6-faq)
7. [Security](#7-security)
8. [License](#8-license)

-----

## 1. Overview

Proteins are large biomolecules that play an important role in the body. Knowing the physical structure of proteins is key to understanding their function. However, it can be difficult and expensive to determine the structure of many proteins experimentally. One alternative is to predict these structures using machine learning algorithms. Several high-profile research teams have released such algorithms, including [OpenFold](https://github.com/aqlaboratory/openfold), [AlphaFold 2](https://deepmind.com/blog/article/alphafold-a-solution-to-a-50-year-old-grand-challenge-in-biology), [RoseTTAFold](https://www.ipd.uw.edu/2021/07/rosettafold-accurate-protein-structure-prediction-accessible-to-all/) ]and others. Their work was important enough for Science magazine to name it the ["2021 Breakthrough of the Year"](https://www.science.org/content/article/breakthrough-2021).

Many AI-driven folding algorithms use a multi-track transformer architecture trained on known protein templates to predict the structure of unknown peptide sequences. These predictions are heavily GPU-dependent and take anywhere from minutes to days to complete. The input features for these predictions include multiple sequence alignment (MSA) data. MSA algorithms are CPU-dependent and can themselves require several hours of processing time.

Running both the MSA and structure prediction steps in the same computing environment can be cost inefficient, because the expensive GPU resources required for the prediction sit unused while the MSA step runs. Instead, using a high-performance computing (HPC) service like [AWS Batch](https://aws.amazon.com/batch/) allows us to run each step as a containerized job with the best fit of CPU, memory, and GPU resources.

This repository includes the CloudFormation template, Jupyter Notebook, and supporting code to run protein analysis algorithms on AWS Batch.

-----

## 2. Quick Start

1. Choose **Launch Stack** and (if prompted) log into your AWS account:

    [![Launch Stack](imgs/LaunchStack.jpg)](https://console.aws.amazon.com/cloudformation/home#/stacks/create/review?templateURL=https://aws-batch-architecture-for-alphafold-public-artifacts.s3.amazonaws.com/main/batch-protein-folding-cfn-packaged.yaml)  
2. For **Stack Name**, enter a value unique to your account and region. Leave the other parameters as their default values and select **Next**.  
    ![Provide a stack name](imgs/name.png)  
3. Select **I acknowledge that AWS CloudFormation might create IAM resources with custom names**.  
4. Choose **Create stack**.  
    ![Choose Create Stack](imgs/create_stack.png)  
5. Wait 20 minutes for AWS CloudFormation to create the necessary infrastructure stack and module containers.  
6. Wait an additional 5 hours for AWS Batch to download the necessary reference data to the Amazon FSx for Lustre file system.  
7. Navigate to [SageMaker](https://console.aws.amazon.com/sagemaker).  
8. Select **Notebook** > **Notebook instances**.  
    ![Select Notebooks](imgs/notebook.png)  
9. Select the **BatchFoldNotebookInstance** instance and then **Actions** > **Open JupyterLab**.
    ![Open BatchFoldNotebookInstance Notebook Instance ](imgs/notebook_instance.png)  
10. Open the quick start notebook at `notebooks/quick-start-protein-folding.ipynb`.  
    ![Open Quick Start Notebook](imgs/open_notebook.png)  
11. Select the **conda_python_3** kernel.  
    ![Open Quick Start Notebook](imgs/select_kernel.png)  
12. Run the notebook cells to create and analyze several protein folding jobs.  
13. (Optional) To delete all provisioned resources from from your account, navigate to [Cloud Formation](https://console.aws.amazon.com/cloudformation), select your stack, and then **Delete**.

-----

## 3. Advanced Configuration

### 3.1. Optional CloudFormation Parameters

- Select "N" for **LaunchSageMakerNotebook** if you do not want to launch a managed sagemaker notebook instance to quickly run the provided Jupyter notebook. This option will avoid the [charges associated with running that notebook instance](https://aws.amazon.com/sagemaker/pricing/).
- Select "N" for **MultiAZ** if you want to limit your Batch jobs to a single availability zone and avoid cross-AZ data transfer charges. Note that this may impact the availability of certain accelerated or other high-demand instance types.
- Provide values for the **VPC**, **Subnet**, and **DefaultSecurityGroup** parameters to use existing network resources. If one or more of those parameters are left empty, CloudFormation will create a new VPC and FSx for Lustre instance for the stack.
- Provide values for the **FileSystemId** and **FileSystemMountName** parameters to use an existing FSx for Lustre file system. If one or more of these parameters are left empty, CloudFormation will create a new file system for the stack.
- Select "Y" for **DownloadFsxData** to automatically populate the FSx for Lustre file system with common sequence databases.
- Select "Y" for **CreateG5ComputeEnvironment** to create an additional job queue with support for G5 family instances. Note that G5 instances are currently not available in all AWS regions.

### 3.2. Manual Data Download

If you set the **DownloadFsxData** parameter to **Y**, CloudFormation will automatically start a series of Batch jobs to populate the FSx for Lustre instance with a number of common sequence databases. If you set this parameter to **N** you will instead need to manually populate the file system. Once the CloudFormation stack is in a CREATE_COMPLETE status, you can begin populating the FSx for Lustre file system with the necessary sequence databases. To do this automatically, open a terminal in your notebooks environment and run the following commands from the **batch-protein-folding** directory:

```python
pip install .
python prep_databases.py
```

It will take around 5 hours to populate the file system, depending on your location. You can track its progress by navigating to the file system in the FSx for Lustre console.

### 3.3. Clean Up

To remove the stack and stop further charges, first slect the root stack from the CloudFormation console and then the **Delete** button. This will remove all resources EXCEPT for the S3 bucket containing job data and the FSx for Lustre backup. You can associate this bucket as a data repository for a future FSx for Lustre file system to quickly repopulate the reference data.

To remove all remaining data, browse to the S3 console and delete the S3 bucket associated with the stack.

-----

## 4. Module Information

### 4.1. JackHMMER

Please visit [https://github.com/EddyRivasLab/hmmer](https://github.com/EddyRivasLab/hmmer) for more information about the JackHMMER algorithm.

### 4.2. AlphaFold

Please visit [https://github.com/deepmind/alphafold](https://github.com/deepmind/alphafold) for more information about the AlphaFold2 algorithm.

The original AlphaFold 2 citation is

```text
@Article{AlphaFold2021,
  author  = {Jumper, John and Evans, Richard and Pritzel, Alexander and Green, Tim and Figurnov, Michael and Ronneberger, Olaf and Tunyasuvunakool, Kathryn and Bates, Russ and {\v{Z}}{\'\i}dek, Augustin and Potapenko, Anna and Bridgland, Alex and Meyer, Clemens and Kohl, Simon A A and Ballard, Andrew J and Cowie, Andrew and Romera-Paredes, Bernardino and Nikolov, Stanislav and Jain, Rishub and Adler, Jonas and Back, Trevor and Petersen, Stig and Reiman, David and Clancy, Ellen and Zielinski, Michal and Steinegger, Martin and Pacholska, Michalina and Berghammer, Tamas and Bodenstein, Sebastian and Silver, David and Vinyals, Oriol and Senior, Andrew W and Kavukcuoglu, Koray and Kohli, Pushmeet and Hassabis, Demis},
  journal = {Nature},
  title   = {Highly accurate protein structure prediction with {AlphaFold}},
  year    = {2021},
  volume  = {596},
  number  = {7873},
  pages   = {583--589},
  doi     = {10.1038/s41586-021-03819-2}
}
```

The AlphaFold-Multimer citation is

```text
@article {AlphaFold-Multimer2021,
  author       = {Evans, Richard and O{\textquoteright}Neill, Michael and Pritzel, Alexander and Antropova, Natasha and Senior, Andrew and Green, Tim and {\v{Z}}{\'\i}dek, Augustin and Bates, Russ and Blackwell, Sam and Yim, Jason and Ronneberger, Olaf and Bodenstein, Sebastian and Zielinski, Michal and Bridgland, Alex and Potapenko, Anna and Cowie, Andrew and Tunyasuvunakool, Kathryn and Jain, Rishub and Clancy, Ellen and Kohli, Pushmeet and Jumper, John and Hassabis, Demis},
  journal      = {bioRxiv}
  title        = {Protein complex prediction with AlphaFold-Multimer},
  year         = {2021},
  elocation-id = {2021.10.04.463034},
  doi          = {10.1101/2021.10.04.463034},
  URL          = {https://www.biorxiv.org/content/early/2021/10/04/2021.10.04.463034},
  eprint       = {https://www.biorxiv.org/content/early/2021/10/04/2021.10.04.463034.full.pdf},
}
```

### 4.3. OpenFold

Please visit [https://github.com/aqlaboratory/openfold](https://github.com/aqlaboratory/openfold) for more information about the OpenFold algorithm.

The OpenFold citation is

```text
@software{Ahdritz_OpenFold_2021,
  author = {Ahdritz, Gustaf and Bouatta, Nazim and Kadyan, Sachin and Xia, Qinghui and Gerecke, William and AlQuraishi, Mohammed},
  doi = {10.5281/zenodo.5709539},
  month = {11},
  title = {{OpenFold}},
  url = {https://github.com/aqlaboratory/openfold},
  year = {2021}
}
```

### 4.4. OmegaFold

Please visit [https://github.com/HeliXonProtein/OmegaFold](https://github.com/HeliXonProtein/OmegaFold) for more information about the OmegaFold algorithm.

The OmegaFold citation is

```text
@article{OmegaFold,
  author = {Wu, Ruidong and Ding, Fan and Wang, Rui and Shen, Rui and Zhang, Xiwen and Luo, Shitong and Su, Chenpeng and Wu, Zuofan and Xie, Qi and Berger, Bonnie and Ma, Jianzhu and Peng, Jian},
  title = {High-resolution de novo structure prediction from primary sequence},
  elocation-id = {2022.07.21.500999},
  year = {2022},
  doi = {10.1101/2022.07.21.500999},
  publisher = {Cold Spring Harbor Laboratory},
  URL = {https://www.biorxiv.org/content/early/2022/07/22/2022.07.21.500999},
  eprint = {https://www.biorxiv.org/content/early/2022/07/22/2022.07.21.500999.full.pdf},
  journal = {bioRxiv}
}
```

### 4.5. RFDesign

Please visit [https://github.com/RosettaCommons/RFDesign](https://github.com/RosettaCommons/RFDesign) for more information about the RFDesign hallucinate and inpainting algorithms.

The RFDesign citation is

```text
@article{RFDesign,
  author = {Jue Wang, Sidney Lisanza, David Juergens, Doug Tischer, Ivan Anishchenko, Minkyung Baek, Joseph L. Watson, Jung Ho Chun, Lukas F. Milles, Justas Dauparas, Marc Expòsit, Wei Yang, Amijai Saragovi, Sergey Ovchinnikov, and David Baker},
  title = {Deep learning methods for designing proteins scaffolding functional sites},
  elocation-id = {2021.11.10.468128},
  year = {2022},
  doi = {10.1101/2021.11.10.468128},
  publisher = {bioRxiv},
  URL = {https://www.biorxiv.org/content/early/2022/07/22/2022.07.21.500999},
  eprint = {https://www.biorxiv.org/content/10.1101/2021.11.10.468128v2.full.pdf},
  journal = {bioRxiv}
}
```

### 4.6. ESMFold

Please visit [https://github.com/facebookresearch/esm](https://github.com/facebookresearch/esm) for more information about the ESMFold algorithm.

The ESMFold citation is

```text
@article{lin2022language,
  title={Language models of protein sequences at the scale of evolution enable accurate structure prediction},
  author={Lin, Zeming and Akin, Halil and Rao, Roshan and Hie, Brian and Zhu, Zhongkai and Lu, Wenting and Smetanin, Nikita and dos Santos Costa, Allan and Fazel-Zarandi, Maryam and Sercu, Tom and Candido, Sal and others},
  journal={bioRxiv},
  year={2022},
  publisher={Cold Spring Harbor Laboratory}
}
```

### 4.7. ProteinMPNN

Please visit [https://github.com/dauparas/ProteinMPNN](https://github.com/dauparas/ProteinMPNN) for more information about the ProteinMPNN algorithm.

The ProteinMPNN citation is

```text
@article{dauparas2022robust,
  title={Robust deep learning--based protein sequence design using ProteinMPNN},
  author={Dauparas, Justas and Anishchenko, Ivan and Bennett, Nathaniel and Bai, Hua and Ragotte, Robert J and Milles, Lukas F and Wicky, Basile IM and Courbet, Alexis and de Haas, Rob J and Bethel, Neville and others},
  journal={Science},
  volume={378},
  number={6615},  
  pages={49--56},
  year={2022},
  publisher={American Association for the Advancement of Science}
}
```

### 4.8. DiffDock

Please visit [https://github.com/gcorso/DiffDock](https://github.com/gcorso/DiffDock) for more information about the DiffDock algorithm.

The DiffDock citation is

```text
@article{corso2023diffdock,
      title={DiffDock: Diffusion Steps, Twists, and Turns for Molecular Docking}, 
      author = {Corso, Gabriele and Stärk, Hannes and Jing, Bowen and Barzilay, Regina and Jaakkola, Tommi},
      journal={International Conference on Learning Representations (ICLR)},
      year={2023}
}
```

### 4.9. RFDiffusion

Please visit [https://github.com/RosettaCommons/RFdiffusion](https://github.com/RosettaCommons/RFdiffusion) for more information about the RFDiffusion algorithm.

The RFDiffusion citation is

```text
@article{joseph_l_watson_broadly_2022,
  title = {Broadly applicable and accurate protein design by integrating structure prediction networks and diffusion generative models},
  url = {http://biorxiv.org/content/early/2022/12/14/2022.12.09.519842.abstract},
  doi = {10.1101/2022.12.09.519842},
  journal = {bioRxiv},
  author = {{Joseph L. Watson} and {David Juergens} and {Nathaniel R. Bennett} and {Brian L. Trippe} and {Jason Yim} and {Helen E. Eisenach} and {Woody Ahern} and {Andrew J. Borst} and {Robert J. Ragotte} and {Lukas F. Milles} and {Basile I. M. Wicky} and {Nikita Hanikel} and {Samuel J. Pellock} and {Alexis Courbet} and {William Sheffler} and {Jue Wang} and {Preetham Venkatesh} and {Isaac Sappington} and {Susana Vázquez Torres} and {Anna Lauko} and {Valentin De Bortoli} and {Emile Mathieu} and {Regina Barzilay} and {Tommi S. Jaakkola} and {Frank DiMaio} and {Minkyung Baek} and {David Baker}},
  year = {2022}
}
```

-----

## 5. Architecture Details

![AWS Batch Architecture for Protein Folding](imgs/batch-protein-folding-arch.png)

### 5.1. Stack Creation Details

This architecture uses a nested CloudFormation template to create various resources in a particular sequence:

1. (Optional) If existing resources are not provided as template parameters, create a VPC, subnets, NAT gateway, elastic IP, routes, and S3 endpoint.
1. (Optional) If existing resources are not provided as template parameters, create a FSx for Lustre file system.
1. Download several container images from a public ECR repository and push them to a new, private repository in your account. Also download a .zip file with the example notebooks and other code into a CodeCommit repository.
1. Create the launch template, compute environments, job queues, and job definitions needed to submit jobs to AWS Batch.
1. (Optional) If requested via a template parameter, create and run a Amazon Lambda-backed custom resource to download several open source proteomic data sets to the FSx Lustre instance.

### 5.2. Cost

There are two types of cost associated with this stack:

- Ongoing charges for data storage, networking, and (optional) SageMaker Notebook Instance usage.
- Per-job charges for EC2 usage and data transfer.

Here are the estimated costs for using the default stack to run [100](https://calculator.aws/#/estimate?id=b1f0310e7266ad45644d3aefaa16f00e11ac1af6) and [5,000](https://calculator.aws/#/estimate?id=8237529417cf8fb38468e6aa8fcf1e6ba9bca527) jobs per month.

To minimize costs, set the `MultiAZ` and `LaunchSageMakerNotebook` options to **N** when creating the stack. This will eliminate the intra-region data transfer costs between FSx for Lustre and EC2 as well as the SageMaker Notebook hosting costs.

-----

## 6. FAQ

Q: When deploying the CloudFormation template, I get an error `Embedded stack arn:aws:cloudformation...  was not successfully created: The following resource(s) failed to create: [AWSServiceRoleForEC2SpotFleetServiceLinkedRole]`. How can I fix this?

This can happen if the service role has already been created in a previous deployment. Try deleting the `AWSServiceRoleForEC2SpotFleetServiceLinkedRole` in the IAM console and redeploy the Cloud Formation template.

-----

## 7. Security

See [CONTRIBUTING](CONTRIBUTING.md#security-issue-notifications) for more information.

-----

## 8. License

This project is licensed under the Apache-2.0 License.

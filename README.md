# PSMA-specific off-the-shelf CAR-NK cells stimulate comprehensive tumoricidal immunity in human metastatic castration-resistant prostate cancer<img src="logo.jpg" width="280px" align="right" />

**Overview**

This repository contains the computational pipeline used for the first-in-human study of healthy-donor-derived, circular mRNA-engineered anti-PSMA CAR-NK cells in patients with metastatic castration-resistant prostate cancer (mCRPC).

The codebase facilitates the integration and analysis of longitudinal, multi-timepoint single-cell transcriptomic and TCR profiling data from 365,183 CAR-NK cells, 2,015,848 peripheral blood cells, and 72,447 tumor-derived cells.

<div align='center' ><b><font size='150'>Overview of CaSee pipeline</font></b></div>
  
<img src="IMG_00001.png" width="3000px"/>


## Pre-requisites:

- Operating System: Ubuntu 22.04.4 LTS

- Central Processing Unit (CPU): AMD EPYC 9754 128-Core Processor

- System Memory (RAM): Approximately 2.1 TB (operating at 5600 MT/s)

- Graphics Processing Unit (GPU): NVIDIA RTX PRO 6000 Blackwell Workstation Edition with 96 GB of VRAM

### Environment and resource allocation

---

For instructions on installing anaconda on your machine (download the distribution that comes with python 3):
https://www.anaconda.com/distribution/

```
~~#conda env create -n CaSee -f configs/CaSee_env_info.yaml~~

# if there are some warnings or errors 
# you can manually install some main packages
# all pip software in all_pip_pacakges.yaml
~~conda create -n CaSee python==3.8.8 # python==3.8~~
conda create -n CaSee python==3.9
conda activate CaSee

pip install pytorch-lightning==1.3.7 # -i https://pypi.tuna.tsinghua.edu.cn/simple
pip install scipy==1.7.0 # -i https://pypi.tuna.tsinghua.edu.cn/simple
pip install numpy==1.20.3 # -i https://pypi.tuna.tsinghua.edu.cn/simple
pip install scanpy==1.7.2 # -i https://pypi.tuna.tsinghua.edu.cn/simple
pip install scikit-learn==0.23.2 # -i https://pypi.tuna.tsinghua.edu.cn/simple
pip3 install opencv-python==4.5.2.54 # -i https://pypi.tuna.tsinghua.edu.cn/simple
pip install torchmetrics==0.3.2 # -i https://pypi.tuna.tsinghua.edu.cn/simple
pip install torchvision==0.10.0 # -i https://pypi.tuna.tsinghua.edu.cn/simple

```

And you also download pytorch https://pytorch.org/ 

Attention, if you in the Chinese mainland, plz use `pip install` instand `conda install` 

**Ubuntu**

```
pip3 install torch==1.9.0+cu111 torchvision==0.10.0+cu111 torchaudio==0.9.0 -f https://download.pytorch.org/whl/torch_stable.html


```

**MacOS**

```
pip3 install torch==1.9.0 torchvision==0.10.0 torchaudio==0.9.0


```

**Windos**

```
pip3 install torch==1.9.0+cu111 torchvision==0.10.0+cu111 torchaudio==0.9.0 -f https://download.pytorch.org/whl/cu111/torch_stable.html
```
> torch==1.9.0+cu111  
> torchvision==0.10.0+cu111  
> torchaudio==0.9.0

# How to install conda/anaconda

This tutorial is for unix/linux users. For the minimal install, we provide tutorial for miniconda installation. Of course, users can download anaconda from [here](https://www.anaconda.com/) and the installation steps are the same.

## Download

### Download through wget

Just copy the following command and paste it in linux shell.

```shell
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

### Download from Browser

Open [miniconda installer page](https://docs.conda.io/en/latest/miniconda.html#linux-installers) and select a version. Here, we download [Miniconda3 Linux 64-bit]() as an example like below.

<center>
    <img style="border-radius: 0.3125em;
    box-shadow: 0 2px 4px 0 rgba(34,36,38,.12),0 2px 10px 0 rgba(34,36,38,.08);" 
    src="./miniconda_download.png">
    <br>
    <div style="color:orange; border-bottom: 1px solid #d9d9d9;
    display: inline-block;
    color: #999;
    padding: 2px;">miniconda download page</div>
</center>

## Install

Use the foloowing command to install miniconda.

```shell
bash Miniconda3-latest-Linux-x86_64.sh
```

Type enter to continue when see this information.

```shell
# just press ENTER
Welcome to Miniconda3 py38_4.9.2

In order to continue the installation process, please review the license
agreement.
Please, press ENTER to continue
>>>
```

Next, you will see User License Agreement, press enetr to finish reading and type "yes".

```shell
Please answer 'yes' or 'no':'
>>> yes
```

Then, decide where to install miniconda. I have a folder named "software" which is already existed in my linux system, and I will install miniconda in this folder.

```shell
Miniconda3 will now be installed into this location:
/home/zhangwei/miniconda3

  - Press ENTER to confirm the location
  - Press CTRL-C to abort the installation
  - Or specify a different location below

[/home/zhangwei/miniconda3] >>> /home/zhangwei/software/miniconda3
```

Here is an important step. I recommened users to initialize Miniconda3 by typing "yes" because we can disable this by modifying .bashrc file. If answering "no", users must enetr miniconda folder to use conda.

```shell
Do you wish the installer to initialize Miniconda3
by running conda init? [yes|no]
[no] >>> yes
```

Now, open a new terminal console, you can see (base) in terminal. Users can install cfDNApipe now.

```shell
# download 
wget https://honchkrow.github.io/cfDNApipe/environment.yml

# clean environment
conda clean -y --all

# install environment
conda env create -n cfDNApipe -f environment.yml

# activate environment
conda activate cfDNApipe

# deactivate environment
conda deactivate

```


## How to disable conda

Assume that users answer "yes" to initialize Miniconda3, when use using a specific environment, some default software will be disabled due to conflicts with the activated virtual environment. There are two ways to disable conda.

### deactivate conda

Users can deactivate conda in (base) environment. Just type conda deactivate in (base) environment.

```shell
(base) [zhangwei@allinone ~]$ conda deactivate
[zhangwei@allinone ~]$
```

### modify bashrc

Modifying bashrc will disable conda permanently. Delete the following lines in ~/.bashrc.

```shell
vi ~/.bashrc
```

Delete or disable the following lines.

``` 
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/zhangwei/software/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/zhangwei/software/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/zhangwei/software/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/zhangwei/software/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
```

Opening a new terminal console, the conda is disabled successfully.

### re-activate conda

Enter the miniconda bin folder.

```shell
cd /home/zhangwei/software/miniconda3/bin
```

Enter the following command.

```shell
[zhangwei@allinone ~]$ cd /home/zhangwei/software/miniconda3/bin
[zhangwei@allinone bin]$ ./conda init
no change     /home/zhangwei/software/miniconda3/condabin/conda
no change     /home/zhangwei/software/miniconda3/bin/conda
no change     /home/zhangwei/software/miniconda3/bin/conda-env
no change     /home/zhangwei/software/miniconda3/bin/activate
no change     /home/zhangwei/software/miniconda3/bin/deactivate
no change     /home/zhangwei/software/miniconda3/etc/profile.d/conda.sh
no change     /home/zhangwei/software/miniconda3/etc/fish/conf.d/conda.fish
no change     /home/zhangwei/software/miniconda3/shell/condabin/Conda.psm1
no change     /home/zhangwei/software/miniconda3/shell/condabin/conda-hook.ps1
no change     /home/zhangwei/software/miniconda3/lib/python3.8/site-packages/xontrib/conda.xsh
no change     /home/zhangwei/software/miniconda3/etc/profile.d/conda.csh
modified      /home/zhangwei/.bashrc

==> For changes to take effect, close and re-open your current shell. <==
```

Opening a new terminal console, the conda is re-activated successfully.

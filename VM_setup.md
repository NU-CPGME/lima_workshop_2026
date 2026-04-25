# Virtual Machine Creation (VirtualBox)

## 1A. Intel/AMD VM

* ISO: lubuntu-24.04.4-desktop-amd64.iso
* Memory: 4096
* Processors: 4
* Storage: 50 GB
* Username: pathogen
* Password: workshop

## 1B: ARM

* ISO: ubuntu-24.04.4-live-server-arm64.iso
* Memory: 4096
* Processors: 4
* Storage: 50 GB
* Username: pathogen
* Password: workshop

> [!IMPORTANT]
> Unselect the "LVM" option when creating storage

```Shell
sudo apt install lxqt cmake
sudo apt clean
reboot
```

## 2: VirtualBox copy/paste setup

* In VirtualBox Host: Settings > General > Features > Shared Clipboard: Bidirectional
* In VirtualBox Client: Devices > Insert Guest Additions CD Image...
* In Lubuntu terminal (Intel/AMD): 

    ```Shell
    sudo /media/pathogen/VBox_GAs_7.2.6/VBoxLinuxAdditions.run
    ```

* In Lubuntu terminal (ARM): 

    ```Shell
    sudo /media/pathogen/VBox_GAs_7.2.6/VBoxLinuxAdditions-arm64.run
    ```

* Reboot

## 3. Install software

```Shell
sudo apt install seaview
```

## 4. Install R

```Shell
sudo apt update -qq
sudo apt install --no-install-recommends software-properties-common dirmngr
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
sudo apt update
sudo apt -y install --no-install-recommends r-base r-base-dev libcurl4-openssl-dev \
    libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
    libfreetype6-dev libtiff5-dev libxml2-dev libcairo2-dev cmake firefox
```

## 5. Install Rstudio 

### Intel/AMD

```Shell
wget https://download1.rstudio.org/electron/jammy/amd64/rstudio-2026.01.1-403-amd64.deb
sudo apt -y install -f ./rstudio-2026.01.1-403-amd64.deb
rm rstudio-2026.01.1-403-amd64.deb
```

### ARM

```Shell
wget https://s3.amazonaws.com/rstudio-ide-build/electron/jammy/arm64/rstudio-2026.04.0-daily-373-arm64.deb  
sudo apt -y install -f ./rstudio-2026.04.0-daily-373-arm64.deb  
rm rstudio-2026.04.0-daily-373-arm64.deb 
```

## 6. Install R packages from terminal

```Shell
sudo Rscript -e 'install.packages(c("tidyverse", "ggsci", "ape", "RColorBrewer", "lubridate","ggpubr", "rlang", "BiocManager")); library(BiocManager); BiocManager::install(c("ggtree", "ggtreeExtra", "treeio"))'
``` 

## 7. Install micromamba

```Shell
"${SHELL}" < <(curl -L micro.mamba.pm/install.sh)
source ~/.bashrc
```

```Shell
micromamba config append channels defaults
micromamba config append channels bioconda
micromamba config append channels conda-forge
micromamba config set channel_priority strict
```

## 8. Create conda environments

### Intel / AMD

```Shell
git clone --depth 1 https://github.com/NU-CPGME/lima_workshop_2026
for env in lima_workshop_2026/conda_envs/*.yaml
do
echo $env
micromamba env create -y -f $env
echo ""
done
micromamba clean -y -a
```

### ARM

```Shell
git clone --depth 1 https://github.com/NU-CPGME/lima_workshop_2026
for env in lima_workshop_2026/conda_envs/ARM/*.yaml
do
echo $env
micromamba env create -y -f $env
echo ""
done
micromamba clean -y -a

mkdir applications
cd applications

## Install snippy 
#### Due to problems with the conda package of perl-bioperl, snippy will not work on ARM machines 
#git clone https://github.com/tseemann/snippy
#echo "export PATH=/home/pathogen/applications/snippy/bin:\$PATH" >> ~/.bashrc

## Install Quast
cd ~/applications
git clone https://github.com/ablab/quast
cd quast
micromamba activate assembly_env
./setup.py install
micromamba deactivate
ln -s quast.py quast
echo "export PATH=/home/pathogen/applications/quast:\$PATH" >> ~/.bashrc

cd ~
source ~/.bashrc
```

## 9. Clean up

```Shell
mv lima_workshop_2026/data/* ~/
cat nanopore_reads/KP.fastq.gz.part-* > nanopore_reads/KP.fastq.gz
rm nanopore_reads/KP.fastq.gz.part-*
if [ -e lima_workshop_2026/example_output ]
then
mv lima_workshop_2026/example_output ~/
fi
rm -rf lima_workshop_2026
```

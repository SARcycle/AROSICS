#!/bin/bash
# -----------------------------------------------------------------------------------
# This script sets up an AWS EC2 instance with AROSICS.
# The instance will be used to create a pre-configured image for use with 
# the "On-Demand Self-Hosted AWS EC2 Runner for GitHub Actions."
#
# Best Practices:
# 1. Ensure the system is fully updated.
# 2. Install all necessary dependencies like Docker, Git, Python, etc.
# 3. Clone and set up AROSICS from the official GitHub repository.
# 4. Configure paths and install additional Python modules.
# 5. The script aims for a clean, minimal setup, leaving the instance ready 
#    for automated tasks triggered by GitHub Actions.
# -----------------------------------------------------------------------------------

# Installing github-runner dependencies
echo "Installing github-runner dependencies..."
sudo apt update -y && \
sudo apt install docker.io -y && \
sudo apt install git -y && \
sudo systemctl enable docker && \
sudo systemctl start docker && \
sudo apt install p7zip-full -y

# Cloning the AROSICS GitHub repository
echo "Cloning GIT..."
git clone https://github.com/SARcycle/AROSICS

# Creating necessary directories
echo "Creating dirs..."
sudo mkdir -p /home/localadmin/
mkdir -p assets
mkdir -p assets/S2

# Installing SARCyle AROSICS installer
echo "Installing SARCyle AROSICS installer..."
sudo bash AROSICS/install_requirements.sh

# Setting paths
echo "Setting paths..."
sed -i 's|/mnt/d/SATROMO/AROSICS_Coregistration/AROSICS/assets/base_data/swissBOUNDARIES3D_1_5_WGS84_buffered_5000m_simplified_DP_5000m.gpkg|/home/ubuntu/base_data/swissBOUNDARIES3D_1_5_WGS84_buffered_5000m_simplified_DP_5000m.gpkg|g' /home/ubuntu/AROSICS/copernicus_api.py
sed -i 's|/mnt/d/SATROMO/AROSICS_Coregistration/AROSICS/assets/base_data/S2_GRI.tif|/home/ubuntu/base_data/S2_GRI.tif|g' /home/ubuntu/AROSICS/coreg_main.py
sed -i 's|/mnt/d/SATROMO/AROSICS_Coregistration/AROSICS/assets/S2|/home/ubuntu/assets/S2|g' /home/ubuntu/AROSICS/util_checkassets.py
sed -i 's|/mnt/d/SATROMO/AROSICS_Coregistration/AROSICS/assets/base_data/log.json|/home/ubuntu/base_data/log.json|g' /home/ubuntu/AROSICS/logger.py
sed -i 's|/mnt/d/SATROMO/AROSICS_Coregistration/AROSICS/assets/base_data/credentials.json|/home/ubuntu/base_data/credentials.json|g' /home/ubuntu/AROSICS/utils.py
sed -i 's|/mnt/d/SATROMO/AROSICS_Coregistration/AROSICS/assets/|/home/ubuntu/assets/|g' /home/ubuntu/AROSICS/main.py

# Installing Python 3.10 venv and modules
echo "Installing python3.10 venv and modules..."
python3.10 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r AROSICS/requirements.txt
pip install earthengine-api
pip install oauth2client
pip install pydrive

# Final instructions
echo "Next step:"
echo "Copy base data from source, (ZIP, encrypt it, transfer to S3 bucket then):"
echo "wget -c <URL to_file>"
echo "7z x -p<password> filename.zip"

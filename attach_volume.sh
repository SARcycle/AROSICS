#!/bin/bash

# Check if volume ID is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <volume-id>"
    exit 1
fi

# Variables
VOLUME_ID=$1  # Get the volume ID from the first argument

# For Older EC2 Instances (Xen-based) /dev/sdX or /dev/xvdX
# For Nitro-based EC2 Instances (e.g., T3, M5, C5, etc.) /dev/sdX
# Nitro-based instances (like T3, M5, C5, etc.): Use /dev/sd[f-p].
# Xen-based instances (like T2, M4, etc.): Use /dev/sd[f-p] or /dev/xvd[f-p].
DEVICE_NAME=/dev/sdf

MOUNT_POINT=/home/ubuntu/assets

#Force detach VOLUME in case it is used by another instance DANGEROUS
# Check the volume state
VOLUME_STATE=$(aws ec2 describe-volumes --volume-ids $VOLUME_ID --query "Volumes[0].State" --output text)

if [ "$VOLUME_STATE" == "in-use" ]; then
    echo "Someone is using the volume: Detaching the volume..."
    aws ec2 detach-volume --volume-id $VOLUME_ID --force

    # Wait for the volume to be detached
    echo "Waiting for the volume to be detached..."
    TIMEOUT=60  # Set timeout in seconds
    elapsed=0   # Initialize elapsed time

    while true; do
        VOLUME_STATE=$(aws ec2 describe-volumes --volume-ids $VOLUME_ID --query "Volumes[0].State" --output>

        if [ "$VOLUME_STATE" == "available" ]; then
            echo "Volume successfully detached."
            break
        fi

        if [ "$elapsed" -ge "$TIMEOUT" ]; then
            echo "Error: Timeout while waiting for volume to detach."
            exit 1
        fi

        sleep 1
        ((elapsed++))
    done
else
    echo "Volume is in state '$VOLUME_STATE'. Detachment not required."
fi




# Get the instance ID
INSTANCE_ID=$(ec2metadata --instance-id)

# Attach the volume
aws ec2 attach-volume --volume-id $VOLUME_ID --instance-id $INSTANCE_ID --device $DEVICE_NAME

# Wait for the volume to be attached
echo "waiting for mount"
TIMEOUT=60  # Set timeout in seconds
elapsed=0   # Initialize elapsed time

while [ ! -e /dev/nvme1n1 ]; do
    if [ "$elapsed" -ge "$TIMEOUT" ]; then
        echo "Error: Timeout while waiting for /dev/nvme1n1 to be available."
        exit 1
    fi
    sleep 1
    ((elapsed++))
done

# Create mount point if it doesn't exist
sudo mkdir -p $MOUNT_POINT

# Mount the volume
sudo mount /dev/nvme1n1 $MOUNT_POINT

echo "done with mount"

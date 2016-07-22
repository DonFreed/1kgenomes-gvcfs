#!/usr/bin/env bash

### Runs on each node upon cluster startup ###
### Configures ephemeral disks             ###

# Install samtools if not installed #
if [[ -z `which samtools` ]]; then
    wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
    tar -jxf samtools-1.3.1.tar.bz2
    cd samtools-1.3.1/
    ./configure
    make install
    cd
fi

# For compute nodes, one of the ephemeral disks is mounted on /mnt #
work=$(df | grep "/mnt" | wc -l)
if [ $work -eq 0 ]; then
    exit 0
fi

# Unmount the alread mounted disk #
mounted=$(mount | grep /mnt | wc -l)
if [ $mounted -eq 1 ]; then
    umount -l /mnt
    # remove /mnt from fstab
    sed -i '/\/mnt/d' /etc/fstab
fi

# Find the ephemeral devices by it's lack of a partition table #
ephemeral=$(fdisk -l 2>&1 1>/dev/null | sed "s/^Disk \([a-Z_:,0-9\/]*\)[a-Z0-9 ,.	']*$/\1/")
n_ephemeral=$(echo $ephemeral | wc -w)

# Create a RAID0 volume #
mdadm --create /dev/md0 --level 0 --raid-devices=$n_ephemeral --run $ephemeral
sleep 5

while [[ `sudo mdadm --detail /dev/md0 | grep 'Rebuild Status'` != '' ]]; do
    sleep 10
done

# Format, mount, set permissions #
mkfs.ext4 /dev/md0
mkdir -p /ephemeral
mount /dev/md0 /ephemeral
chmod ugo+rwx /ephemeral

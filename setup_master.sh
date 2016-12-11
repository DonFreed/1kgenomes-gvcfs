#!/usr/bin/env bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cp ${DIR}/sample_fastq_to_gvcf.py /data/
echo 'PATH=/opt/sge6/bin/linux-x64/:$PATH' >> /home/onekg/.bashrc

# Configure ephemeral disk space and RAM as a SGE resource (complex) #
qconf -sc | head -n -1 > ~/tmp_complex.conf
echo "ephemeral           ephem      MEMORY      <=    YES         YES        0        0" >> ~/tmp_complex.conf
echo "total_mem           t_mem      MEMORY      <=    YES         YES        0        0" >> ~/tmp_complex.conf
qconf -Mc ~/tmp_complex.conf

# What are the cluster hosts? #
hosts=$(qconf -sel)
n_hosts=$(echo $hosts | wc -w)

# Run the ephemeral configuration script on each node #
for host in $hosts
do
    ssh $host "bash -s" < ${DIR}/node_startup.sh &
done

# Mount ephemeral disks on the compute nodes #
nodes_up=0
while [ $nodes_up -lt $(($n_hosts - 1)) ]
do
    sleep 10
    nodes_up=0
    for host in $hosts
    do
        is_up=$(ssh $host df | grep "ephemeral" | wc -l)
        nodes_up=$((nodes_up + is_up))
    done
done

# Get the ephemeral space for each node #
for host in $hosts
do
    # Get the inital configurations #
    qconf -se $host | grep -v "complex_values\|load_values\|processors\|^ " > ~/tmp_node.conf

    # Find the ephemeral space (in bytes) #
    host_size=$(ssh $host df -B 1 | grep "ephemeral" | awk '{print $2}')

    # Find the total memory (in bytes) #
    mem_size=$(ssh $host cat /proc/meminfo | grep "MemTotal" | sed 's/^MemTotal:[ ]*\([0-9]*\) kB$/\1/')
    mem_size=$((mem_size * 1024))

    if [ -z $host_size ]; then
        echo "complex_values        ephemeral=0,total_mem=$mem_size" >> ~/tmp_node.conf
    else
        echo "complex_values        ephemeral=$host_size,total_mem=$mem_size" >> ~/tmp_node.conf
    fi

    qconf -Me ~/tmp_node.conf
done

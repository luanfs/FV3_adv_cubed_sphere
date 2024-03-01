#!/bin/bash
# Luan Santos - 2023
# Script to get data from remote ybytu servers.

#-------------------------------------------------------------------------------------------------------
# Hosts and users
remote_host1="brucutu.ime.usp.br"            # accessed from local host
remote_host2="ybytu"      # accessed from remote host1
user1="luansantos"                   # user at remote_host1
user2="luansantos"                   # user at remote_host2
codedir="FV3_adv_cubed_sphere"
#-------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
# Directories and output
output_dir_remote_host2="/home/"$user1"/$codedir/" #remote host data directory (remote_host2)
output_dir_local_host="/home/luanfs/doc/$codedir/"  #local data directory
data='data'           #name of directory where data is in output_dir_remote_host2
output=$data   #output file (.tar)

#-------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
# Fisrt, we create a tarball with the target data in remote host 2
#-------------------------------------------------------------------------------------------------------
echo "Creating tarball at $remote_host2:"
ssh -t -CYttg -L 10000:localhost:10000 $user1@$remote_host1 ssh -CYtt -L 10000:localhost:22 $user2@$remote_host2<<EOF
	cd $output_dir_remote_host2; pwd; ls; tar -cvf $data.tar $data; exit
EOF
echo "Created tarball at $remote_host2"

#-------------------------------------------------------------------------------------------------------
# Get data from remote_host1 to local
echo
echo "Getting data from $remote_host2 to local:"
echo $output_dir_remote_host2$output 
echo $output_dir_local_host
 rsync -e "ssh -o ProxyCommand='ssh -A -W %h:%p $user1@$remote_host1 ssh -A -W %h:%p $user2@$remote_host2'" -avz luansantos@ybytu:$output_dir_remote_host2/$data $output_dir_local_host/
echo "Got data from $remote_host2 to local"
echo
#-------------------------------------------------------------------------------------------------------


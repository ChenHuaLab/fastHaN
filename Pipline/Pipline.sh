#!/bin/bash

Usage()
{ echo -e "NAME"
    echo -e "\t$0 (from unaligned fasta to haplotype network)"
    echo -e "USAGE"
    echo -e "\tbash $0 -i <input.fasta> -o <output_prefix> -n <network_type[mjn|msn|modified_tcs|original_tcs]> -t <thread_number>"
    echo -e "CAUTION"
    echo -e "\tEach of the parmeters above should be provided!\n"
}


Option()
{
    while [ -n "$1" ]; do
    case "$1" in
        -i)
            In=$2
            echo -e "\t[Option] Found the path of input fasta file ($2)"
            ;;
        -o)
            Out=$2
            echo -e "\t[Option] Found the prefix of output graph file ($2)"
            ;;
        -n)
            Type=$2
            echo -e "\t[Option] Found the type of haplotype netowrk ($2)"
            ;;
        -t)
            Threads=$2
            echo -e "\t[Option] Found the thread number used in network construction ($2)"
            ;;
         *)
            echo -e "\t[Option] $1 is not a valid option, please check the usage\n"
            Usage
            exit -1
            ;;
    esac
    shift 2
    done
}


Main()
{
    if [ $# -lt 8 ];then
        Usage
        exit -1
    else
        echo -e "[`date`] parse the parameters of command line ..."
        Option $*
        echo ""
    fi

    # generate the alignment by muscle3 (must be executable)
    echo "[`date`] generate the alignment by muscle3 ..."
    muscle3 -in ${In} -out ${Out}".alignment"

    # convert the alignment to phylip file by Fasta2Phylip.py (must be executable)
    echo "[`date`] convert the alignment to phylip file by Fasta2Phylip.py ..."
    python Fasta2Phylip.py ${Out}".alignment" > ${Out}".phylip"

    # construct the haplotype network by fastHaN (must be executable)
    echo "[`date`] construct the haplotype network by fastHaN ..."
    case "${Type}" in
        "mjn")
            ./fastHaN mjn -i ${Out}".phylip" -t ${Threads} -o ${Out}
            ;;
        "msn")
            ./fastHaN msn -i ${Out}".phylip" -o ${Out}
            ;;
        "modified_tcs")
            ./fastHaN modified_tcs -i ${Out}".phylip" -t ${Threads} -o ${Out}
            ;;
        "original_tcs")
            ./fastHaN original_tcs -i ${Out}".phylip" -t ${Threads} -o ${Out}
            ;;
         *)
            echo -e "[fastHaN] $1 is not a valid network type, please check the usage\n"
            Usage
            exit -1
            ;;
    esac

    echo "[`date`] done!"

}

# entrance of the pipline
Main $*

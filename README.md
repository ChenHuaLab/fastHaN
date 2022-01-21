fastHaN
=========================
A fast and scalable program for constructing haplotype network for large-sample sequence data sets. <br>
You can also get the source code of fastHaN at http://chenlab.big.ac.cn/software/submitInfo.html.


Description
=========================
* fastHaN is a fast and scalable program for constructing haplotype networks for large samples <br>
* fastHaN can implement the minimum joint network (MJN) and the Templeton-Crandall-Sing (TCS) algorithms <br>
* The implementation in single-thread mode is much faster than the existing softwares <br>
* Furthermore, fastHaN enables multi-threaded mode with good scalability


Quick start guide
========================
    For Linux:
    1. download the software and grant execution permissions
       chmod +x fastHaN_linux

    2. run the test data (Example/Test1000.phy.gz) with mjn algorithm
       ./fastHaN_linux mjn -i Example/Test1000.phy.gz -t 8 -o Test1000.gml

    For windows:
    1. run in the CMD window
        fastHaN_win.exe mjn -i Example/Test1000.phy.gz -t 8 -o Test1000.gml


Usage
========================
    fastHaN [Options]
    
    Options:
        original_tcs        original TCS algorithm (Clement et al. 2002)
        modified_tcs        optimization of TCS implementated by PopART (Leigh and Bryant, 2015)
        msn                 optimization of MSN implementated by PopART
        mjn                 optimization of MJN implementated by PopART
         

Options
========================

### original_tcs

    usage: fastHaN original_tcs [arguments]

    input:
        -i    input phylip format file
    options:
        -t    (int)thread number(default:8)
        -a    (int)is mark the sites which contains ambiguous base, 1:mask, 0:not(default:0)
        -m    (int)is merge intermediate vertex, 1:merge, 0:not:(default:0)
    output:
        -o    output graph file (*.gml)
            

### modified_tcs

    usage: fastHaN modified_tcs [arguments]
    
    input:
        -i    input phylip format file
    options:
        -t    (int)thread number(default:8)
    output:
        -o    output graph file (*.gml)
            
            
### msn

    usage: fastHaN msn [arguments]
         
    input:
        -i    input phylip format file
    options:
        -e    (int)epsilon(default:0)
    output:
        -o    output graph file (*.gml)
            
            
### mjn

    usage: fastHaN mjn [arguments]
         
    input:
        -i    input phylip format file
    options:
        -t    (int)thread number(default:8)
        -e    (int)epsilon(default:0)
    output:
        -o    output graph file (*.gml)

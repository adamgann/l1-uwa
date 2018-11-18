This package consists of Matlab codes set_channel_params.m and channel_simulator.m and 
a series of Matlab functions that will be called by the main codes. 

To change the channel parameters, edit set_channel_params.m and compile to generate the
[file_name].prm and [file_name].dop data files. Alternatively, you can directly 
update the data files via a text editing software. 
Once a set of data files are generated, run channel_simulator.m to simulate the channel 
for [file_name].prm and [file_name].dop.
 
The simulator takes into account both large- and small-scale variations as well as motion-
induced Doppler. For more information on the channel model, 
please visit http://millitsa.coe.neu.edu/?q=research/simulator

For more information about the Takagi factorization method and the enclosed implementation
package, please visit http://www.cas.mcmaster.ca/~qiao/software/takagi/matlab/readme.html


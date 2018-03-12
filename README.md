## Introduction - E-Fi overview
This project contains the base code for the E-Fi proposal, an ongoing research from Northeastern University in the Genesys Lab in collaboration with the Mathworks inc. E-Fi allows Wi-Fi devices to coexist with LTE on the same band without requiring any changes on the LTE and Wi-Fi Physical Layer (unlike Listen-Before-Talk or CSAT) by extending its capabilities with the WiFi Direct Standard, which allows for a multi-hop communication, and scheduling its transmissions within Almost Blank Subframes (ABS) (an Inter-cell Interference coordination mechanism present in LTE Release 10). The E-Fi mechanism can be summarized into the following steps:

1. Packet Success Rate (PSR) evaluation based on received power (Prx) and Signal-to-Noise Ratio (SINR).
2. Node Pre-categorization (Relay, Wi-Fi Direct Client or Regular Client).
3. Device discovery and PSR exchange.
4. Node Grouping and ABS distribution.

## Project hierarchy

The root folder contains several subfolders that can be categorized as:

1. **Data** (*/0_DATA*, serves as a pre-requisite of the E-Fi experiments)
2. **Simulation** (*/1_TABLE_PHY* and */2_TABLE_PHY*). They generate the tables needed in E-Fi experiments. Sample tables are already available in /0_DATA.
3. **E-Fi experiments** (the rest). These experiments require the data stored under */0_DATA*.

In addition, */Configuration_MAC.m* and */Configuration_PHY.m* are the configuration files for the *simulations* and *E-fi experiments*.  For experiments 4 (/4_EXPERIMENT_NODE_CATEGORY), 5 (/5_EXPERIMENT_FEASIBLE_REGION) and 6 (/6_EXPERIMENT_GENERAL_STATS), the code calls the function *EFi_function.m*:

- **EFi_function.m**: This function emulates E-fi as its best. It randomly deploys Wi-Fi nodes within the coverage area, calculate their Receive Power and SINR based on its location and the location of the BS and AP, calculates the PSR based on the TABLE_PER_ABS0/1, calculates the PSR between the Relay candidates and Wi-Fi Direct candidates and group the nodes based on the output by the Hungarian algorithm. The code returns the node categorization, the location vector and measurements of the Initial/final PSR and Throughput based on the TABLE_NEW_MAC.

**NOTE**: *For the correct execution of the code, ensure you are running it on version R2016a or later with the LTE System Toolbox, WLAN System toolbox and Parallel Toolbox. The later one can be deactivated in the code, though it would increase the execution time.*

### 1. Create PHY Table
The code is under */1_TABLE_PHY* directory and uses the LTE and WLAN System toolboxes to perform the mapping PSR - <Prx,SINR> (PSR_TABLE). This information is supposed to be known by every Wi-Fi device once E-fi starts operating. The output of the code is a 3D Table, where each combination of received power and SINR map to a PSR value. The code to characterize the PSR is:

- **LTE_WiFi_TXChRx_Runnable.m**: This is the main script. It contains the LTE and Wi-Fi configurations (tunable) and generates the testing grid by iterating over different values of AP- BS distances and node location. For each of this values, the script calls *LTE_WiFi_TXChRx.m*, which returns metrics on the Wi-Fi and LTE performance.
- **LTE_WiFi_TXChRx.m**: This code evaluates the throughput, BER and Frame Detection Errors for LTE and Wi-Fi standard when coexisting using Almost Blank Subframes (ABS). Wi-Fi transmissions are always scheduled in ABS, in which LTE only transmits Control Signal. LTE Data is scheduled in non-ABS. The performance is evaluated only at a specified location (distance from the BS and AP).

### 2. Create MAC Table
The code is under */2_TABLE_MAC* directory and takes [1] as a baseline to model the MAC procedure using a state machine. The input of the function is the PSR_TABLE (generated in */1_TABLE_PHY*) and it provides the probability of dropping a packet, the number of attempts to transmit a packet and the average time to transmit. The later, combined with the length of the PSDU gives us an estimation of the throughput. The output is defined as TABLE_MAC and is used in most of the further experiments to compute the Throughput accounting for retrasnmissions following the DCF (Distributed Cordination Function) present in the 802.11 standards family.

- **CreateTableMAC.m**: Main executable. The code runs for several values of *PSDU Length* and *alpha*. *Alpha* is a parameter that modifies the contention window in E-Fi to allow for an evenly distributed channel access between AP's and Relays (e.g. Alpha equals to 0 gives Relays more chances to access the channel than the AP to forward packets to Wi-Fi Direct Clients). *PSDU Length* represents the length in bits of the PSDU field in a packet. 

### 3. Throughput Experiment
The simulation set is under */3_EXPERIMENT_THROUGHPUT* directory and generates part of the TABLE_NEW_MAC and shows the results. Although not necessary for E-Fi performance evaluation, it serves as a good visualization of the table creation mechanism.

- **ThroughputExp.m**: Main executable. This code is a particularization of **CreateTableMAC.m** (in /2_TABLE_MAC). The code evaluates the Throughput for a range of values of either *PSDU Length* or *alpha*, keeping the other parameter to a fixed value throughout the simulation.
≤
### 4. Node Categorization Experiment
The simulation set is under */4_EXPERIMENT_NODE_CATEGORY* and helps to understand how nodes are statistically categorized as Group Owners/Relays (GO), Wi-Fi Direct Clients (WDC) or Clients (CSZ and CNSZ) in the X-Y plane (spatial categorization based on the total number of available nodes in the area). 

- **NodeCategorization.m**: Main executable. The Wi-Fi nodes are located randomly within the coverage area. The code runs over a large number of iterations to ensure we cover most locations and we obtain meaningful results (200 or more iterations are recommended). The output shows (1) the node characterization and (2) the PSR improvement for the Wi-Fi Direct Clients by employing E-Fi.

### 5. Feasible Region Experiment
The simulation set is under */5_EXPERIMENT_FEASIBLE_REGION* and shows how the quality of the communications between the AP and the Wi-Fi Direct Client (PSR between Ap and WDC), the GO and the WDC (PSR between GO and WDC) and the AP and the GO (PSR between AP and GO) impact on the candidacy analysis and final node grouping in E-Fi. 

- **FeasibleRegion.m**: Main executable. The code will run for several iterations, select the combinations (PSR between Ap and WDC) that meet the initial number of transmission (within certain margins of tolerance) and compare them with the ones E-Fi provides. The 2D output helps to understand the distribution and final categorization. The 3D output also show the estimation of the Throughput for the selected PSR configurations. To that end, a list of initial number of transmissions for success needs to be defined. This can be seen as a list of 1/(PSR between Ap and WDC), since this is the initial quality of the link without E-Fi.

### 6. Global Experiment
The simulation set is under */6_EXPERIMENT_GENERAL_STATS* and generates an extensive analysis on the performance of E-Fi in terms of PSR and Throughput as a function of the Number of Nodes (Nnodes), the PSR Threshold (PERth, defining the Safe Zone in E-Fi, or PSRth, being PSRth = 1 - PERth), the BS location (BSloc, distance between AP and BS) and Alpha (modifies the contention window in E-Fi). 

- **GeneralStats.m**: Main executable. The parameter “exp” within the executable *GeneralStats.m* controls the parameter to evaluate in the experiment and can be set to: 'Nnodes', 'PERth', 'PSRth', 'BSloc' or 'Alpha’. This is the most important experiment set in E-Fi. Once the execution finalizes, it stores the results in *dataGeneralStats.m* and calls GeneralStatsPlot.m.
- **GeneralStatsPlot.m**: Plots the results stored in the mat file *dataGeneralStats.mat*. It prints a timestamp from the previous execution as a reference.
- **dataGeneralStats.mat**: Data generated after executing *GeneralStats.m*. It contains the relevant parameters to evaluate the performance and plot the results. It changes everytime the script *GeneralStats.m* is executed.

## Contact

For more information about the code or E-Fi, please don't hesitate to contact us by email. We'll be happy to assist you and answer any question you may have.

Carlos Bocanegra Guerra  
PhD Candidate  
Electrical and Computer Engineering (EECE)  
Northeastern University  
bocanegrac@coe.neu.edu

Zhengnan Li  
Master of Science candidate  
Electrical and Computer Engineering (EECE)  
Northeastern University  
li.zhengn@husky.neu.edu  

Kaushik R. Chowdhury  
Associate Professor  
Electrical and Computer Engineering (EECE)  
Northeastern University  
krc@ece.neu.edu

## Acknowledgements

This work is supported in part by MathWorks under the Development-Collaboration Research Grant and by the U.S. Office of Naval Research under grant number N00014-16-1-2651. We would like to thank Mike McLernon, Ethem Sozer, Rameez Ahmed and Kunal Sankhe for their continued support and guidance on this project.

## References

[1] R. Subramanian, B. Drozdenko, E. Doyle, R. Ahmed, M. Leeser and K. R. Chowdhury, "High-Level System Design of IEEE 802.11b Standard-Compliant Link Layer for MATLAB-Based SDR," in IEEE Access, vol. 4, no. , pp. 1494-1509, 2016.
Code:  https://github.com/80211bSDR. 
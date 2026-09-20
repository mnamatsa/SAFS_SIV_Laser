# DOUBLE PULSE ONLINE DMD
# Description
This is a online DMD code based on the Zhang, H., Rowley, C.W., Deem, E.A., Cattafesta, L.N.: Online dynamic mode decomposition for time-varying systems. SIAM Journal on Applied Dynamical Systems 18(3), 1586–1609 (2019) https://doi.org/10.1137/18M1192329
It was adapted to work with double pulse system and can resolve Dynamic Modes up to frequency 1/2dt with dt being delay between two double pulses. It currently utilizes a Phantom MATLAB SDK to process .cine files directly
# Phanto SDK
The SDK can be requested here: https://phantomhighspeed.my.site.com/PhantomCommunity/s/article/Phantom-SDK-MATLAB-LABVIEW-Access
# Results
 M=0.98 jet (dt=800 ns) | Boundary Layer of a M=0.91 plug nozzle (dt=800 ns) |
|---|---|
| ![Shear Layer DMD Mode f = 74kHz Hz](https://raw.githubusercontent.com/mnamatsa/SAFS_SIV_Laser/blob/main/OnlineDMD/ODMD%20f%3D74349.8%20Hz.gif?raw=true) | ![Boundary Layer DMD Mode f= 315Khz](https://raw.githubusercontent.com/mnamatsa/SAFS_SIV_Laser/blob/main/OnlineDMD/ODMD%20f%3D315097.8%20Hz.gif?raw=true) |

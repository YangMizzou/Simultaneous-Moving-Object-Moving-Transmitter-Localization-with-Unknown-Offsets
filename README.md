# Simultaneous-Moving-Object-Moving-Transmitter-Localization-with-Unknown-Offsets
MATLAB Processing Code for Multistatic Localization of a Moving Object by a Moving Transmitter at Unknown Location Having Unknown Time and Frequency Offsets

If you use any of the following codes in your research, please cite the corresponding paper as a reference in your publication. Thank you!

# Project Abstract

This project investigates the multistatic localization of a moving object in position and velocity using an uncoordinated moving transmitter of unknown position and velocity. The measurements are time delays and frequencies and each kind is subject to an unknown amount of offset. We have shown that incorporating the direct-path measurements between the transmitter and receivers to the indirect-path measurements from the transmitter through the object to the receivers improves the localization performance, although nuisance parameters including the transmitter position, velocity and offsets are required for estimation. The condition about the localization geometry that can eliminate the degradation due to the offsets is derived for IID Gaussian noise. Algebraic solution to localize the moving object is developed together with the performance analysis in reaching the Cramer-Rao Lower Bound accuracy under Gaussian noise over the small error region. The particular case of having time delay measurements only is examined and the optimal geometry for handling unknown transmitter position and time offset is devised. Simulations validate the theoretical developments.

# Code Description

Multistatic Localization of a stationary object by a stationary transmitter at unknown position having unknown time offset:

MSLocJntObjTxPos.m: Algebraic Closed-Form Solution
MSLocJntObjTxPos_MLE.m: Iterative MLE Solution 
MSLocObjPosInd_MLE.m: Iterative MLE Solution by Indirect-path measurements
MSLocObjPosCRLB.m: CRLBs
Example_Fig6Fig8.m: Fig 6 and 8 Reproduction

Multistatic Localization of a moving object by a moving transmitter at unknown location having unknown time and frequency offsets:

MSLocJntObjTxPosVel.m: Algebraic Closed-Form Solution
MSLocJntObjTxPosVel_MLE.m: Iterative MLE Solution
MSLocObjPosVelInd_MLE.m: Iterative MLE Solution by Indirect-path measurements 
MSLocObjPosVelCRLB.m: CRLBs
Example_Fig4Fig7: Example

# Reference

 Y. Zhang and K. C. Ho, "Multistatic moving object localization by a moving transmitter of unknown location and offset," IEEE Trans. Signal Process., vol. 68, pp. 4438-4453, 2020.
 
 https://cisp.ece.missouri.edu/index.html

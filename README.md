## Summary
This code is based on the paper 
A. Shahmansoori, G. E. Garcia, G. Destino, G. Seco-Granados and H. Wymeersch, "Position and Orientation Estimation Through Millimeter-Wave MIMO in 5G Systems," in *IEEE Transactions on Wireless Communications*, vol. 17, no. 3, pp. 1822-1835, March 2018.

The matlab code (main.m) generates a 2D environment with random scatterers and line-of-sight (LOS) between a transmittee (at a fixed location) and a receiver (with unknown location and orientation). Both the transmitter and receiver are equipped with uniform linear arrays (for the transmitter the ULA is aligned with the vertical axis). The transmitter sends a sequence of beams, with associated precoders, to the receiver. At the receiver, the time-of-arrival (TOA), angle-of-arrival (AOA), and angle-of-departure (AOD) are estimated using distributed compressed sensing (DCSSOMP.m) in the beamspace domain. The extended version (extended_main.m) is an improvement of main.m with a refinement step using Golden-Section search for AOA, AOD and Least Square (LS) for channel amplitude and TOA.

## Authors

The code was developed by:
* **[Henk Wymeersch](https://sites.google.com/site/hwymeers/)**
* **[Gonzalo Seco Granados](http://spcomnav.uab.es/~gseco/)**

The extended version (extened_main.m) with a refinement step was developed by: 
Ngoc-Son Duong (sondn24@vnu.edu.vn), Quoc-Tuan Nguyen, Thai-Mai Dinh-Thi (UET-VNU)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

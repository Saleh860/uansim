# UANSim

The Underwater Acoustic Network Simulator (UANSim) is a simulator of underwater network protocols designed to allow researchers with MATLAB programming experience to experiment with existing underwater network protocols.

## 1. Supported Protocol Stack

1. Application Layer Protocols
   - Constant bitrate data source
2. Network Layer Protocols
   - Flooding
   - Depth-based routing (DBR)
   - Resilient pressure-based routing (RPR)
   - Depth-based probabilistic routing (DPR)
3. Link Layer Protocols
   - ALOHA
4. Physical Layer Protocols
5. Physical Channel Model
   - Thorp's propagation model (single path, spreading loss, frequency-dependent absorption loss)
   - Ideal propagation model (single path, no path loss)
   - Ambient channel noise (turbulance, shipping, surface waves, and thermal)
   - Additive white Gaussian noise model
   - Support multipath propagation using Power Delay Profiles
     The simulation configuration is controlled through a MATLAB script.

// This file contains the documantation description of some of the functions used in getting the truth IDE (energy deposits) from simualtion

// This is a vector which conatins information about each channel.
std::vector<const sim::SimChannel*> simChannelVec_t 

event.getView(fSimulationProducerLabel, simChannelVec); // ++ Find out what this line does. 

// A vector of particles: This object keeps track of the relevant information for a particle that is to be tracked through the simulation. 
// It has data members for position and momentum along a particle's trajectory as well as the mother and daughter information.
art::Handle<std::vector<simb::MCParticle>> particleHandle; 

/* Sim::SimChannel
Energy deposited on a readout channel by simulated tracks.
This class stores the list of all energies deposited on a readout channel. The number of electrons is stored as well.
The information is organized by time: it is divided by TDC ticks, and each TDC tick where some energy was deposited appears in a separate entry, while the quiet TDC ticks are omitted. For each TDC, the information is stored as a list of energy deposits; each deposit comes from a single Geant4 track and stores the location where the ionization happened according to the simulation (see sim::IDE class).
Note that there can be multiple energy deposit records (that is sim::IDE) for a single track in a single TDC tick.
*/
sim::SimChannel


//Type of list of energy deposits for each TDC with signal.
sim::SimChannel::TDCIDEs_t

/*
Ionization at a point of the TPC sensitive volume.
This class stores information about the ionization from the simulation of a small step of a track through the TPC active volume.
Ionization information consists of both energy and number of electrons. It is of paramount importance to understand what each field stores:

position: where the ionization occurred (from Geant4 simulation)
track ID: Geant4 track ID of the ionizing particle
energy: amount of energy released by ionization (from Geant4 simulation)
electrons: amount of electrons reaching the readout channel
Note the different definition of the electrons respect to the rest: it describes the electrons at the anode after the drifting occurred, while all the other quantities can be related to the moment the ionization happened.

The number of electrons typically includes inefficiencies and physics effects that reduce and spread the electrons. In the simulation, this yields a fractional number of electrons.
Each IDE is also typically associated with a time (TDC) count, that is the time at which the ionized electrons reached the readout channel, in electronic ticks, as opposed as the time when ionization occurred. The latter is not stored.
At the time of writing this documentation (LArSoft 6.4.0), IDEs are computed in larg4::LArVoxelReadout. The energy and track ID come directly from Geant4 simulation. The position is the mid point of the Geant4 step that produced ionization. The electrons are

converted from that same energy (using a fundamental conversion factor stored in larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h)
applied recombination effect by larg4::IonizationAndScintillation::Reset()
applied attenuation and diffusion in larg4::LArVoxelReadout::DriftIonizationElectrons()
The latter also assembles the sim::IDE objects to be stored into sim::SimChannel
*/
sim::IDE

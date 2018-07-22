
typedef std::vector<const sim::SimChannel*> simChannelVec_t; // vector of channels


simChannelVec_t simChannelVec; // Decalre simChannleVec as a simChannelVec type
event.getView(fSimulationProducerLabel, simChannelVec);

// create vector of particles involved in event
art::Handle<std::vector<simb::MCParticle>> particleHandle;
event.getByLabel(fSimulationProducerLabel, particleHandle);

double dEnergyInitial = 0.;
// find primary particle to get initial truth energy
for (std::vector<simb::MCParticle>::const_iterator particlePtr = particleHandle->begin(); particlePtr != particleHandle->end(); ++particlePtr)
{
    const simb::MCParticle &particle = (*particlePtr); // De-reference the particle pointer
    if (particle.Process() == "primary")
    {                                                       // if particle is primary
        dEnergyInitial += particle.Momentum(0).E() * 1000.; // initial energy (adds together over all primary particles)
    }
}

double dEnergyDepositedSum = 0.; // total amount of deposited energy

// loop over channels
for (simChannelVec_t::const_iterator channelPtr = simChannelVec.begin(); channelPtr != simChannelVec.end(); ++channelPtr)
{
    // get individual channel
    const sim::SimChannel &channel = *(*channelPtr);

    // select plane 0
    if (fGeometry->SignalType(channel.Channel()) == geo::kCollection)
    {
        // get time slice
        const sim::SimChannel::TDCIDEs_t &timeSlices = channel.TDCIDEMap();

        // loop over time slices which contain a signal
        for (sim::SimChannel::TDCIDEs_t::const_iterator timePtr = timeSlices.begin(); timePtr != timeSlices.end(); ++timePtr)
        {

            // get vector of IDEs (second entry in map<first,second>)
            typedef std::vector<sim::IDE> depositVec_t;
            const depositVec_t &energyDeposits = (*timePtr).second;

            // loop over IDEs
            for (depositVec_t::const_iterator energyPtr = energyDeposits.begin(); energyPtr != energyDeposits.end(); ++energyPtr)
            {
                // get IDE
                const sim::IDE &energyDeposit = (*energyPtr);
                TruthEnergyDeposits.push_back(energyDeposit.energy);
                TruthZPos.push_back(energyDeposit.z - 250.);  // add z pos to a tree
                TruthXPos.push_back(energyDeposit.x - 102.5); // add xpos to tree

                // Add hit coordinates to PCA vectors

                if (use3D == true)
                {
                    TVectorD vHitPos(3);
                    vHitPos[0] = energyDeposit.x - 102.5; // X
                    vHitPos[1] = energyDeposit.y;         // Y
                    vHitPos[2] = energyDeposit.z - 250.;  // Z

                    if (vHitPos[2] < 0.)
                    {
                        vHitPos[0] = 0.;
                        vHitPos[1] = 0.;
                        vHitPos[2] = 0.;
                        addDataEntry(0, vHitPos);
                        dEnergyDepositedSum += 0;
                    } // delete any spurious data points (z<0)
                    else
                    {
                        addDataEntry(energyDeposit.energy, vHitPos);
                        dEnergyDepositedSum += energyDeposit.energy;
                    } // Add the entries to the vector

                    //addDataEntry(energyDeposit.energy, vHitPos); // Add the entries to the vector
                    //dEnergyDepositedSum += energyDeposit.energy; // calculate energy sum
                    hTruePosition3D->Fill(vHitPos[0], vHitPos[1], vHitPos[2], energyDeposit.energy); // fill 3D PCA histogram
                }
                else
                { // 2D case

                    TVectorD vHitPos(2);
                    vHitPos[0] = energyDeposit.x - 102.5; // X
                    vHitPos[1] = energyDeposit.z - 250.;  // Z

                    if (vHitPos[1] < 0.)
                    {
                        vHitPos[0] = 0.;
                        vHitPos[1] = 0.;
                        addDataEntry(0, vHitPos);
                        dEnergyDepositedSum += 0;
                    } // delete any spurious data points (z<0)
                    else
                    {
                        addDataEntry(energyDeposit.energy, vHitPos);
                        dEnergyDepositedSum += energyDeposit.energy;
                    } // Add the entries to the vector

                    //dEnergyDepositedSum += energyDeposit.energy; // calculate energy sum

                    if (HistFill != event)
                    {
                        hTruePosition2D->Fill(vHitPos[1], vHitPos[0], energyDeposit.energy);
                    }
                }

            } // End loop over energy deposits (IDEs)

        } // End loop over Time slices

    } // End condition if in collection plane

} // End loop over all channels

# Start/Set Parameters
  showerstart = 100     #Punch-Through Start Assumption
  shower = false  #Boolean Shower Condition
  mipLimitBeam = 6.0 + 0.1*beamEnergy     # Energy Threshold
  hitLimitBeam = int(3.77 + 1.44*log(beamEnergy) + 0.5)   #Hit Treshold
  window = 6
  first_layer = 1
  last_layer = 39


  mip0 = mip1 = 1.4
  hit0 = hit1 = 1     # Values for virtual MIP layers before first layer of calorimeter to calculate window average, assumed MIP track with 1 hit and energy 1.4 MIP

  i = first_layer - 1  // begin search

  while (  i < last_layer  and  !shower ):

    #Virtual Layer Case
    if  i < window:
        mip1 = mip0 + ( energy_in_layers[i+1] - 1.4 )/window
        hit1 = hits_in_layers[i+1]


    else:                            #Normal Layer Case
      mip1 = mip0 + ( energy_in_layers[i+1] - energy_in_layers[i+1-window])/window;
      hit1 =hits_in_layer[i+1]


    if ( (mip0 + mip1) > mipLimitBeam  and  ( hit0 + hit1 ) > hitLimitBeam ):
       shower = true


     mip0 = mip1
     hit0 = hit1
     i++

 if shower == true:
    showerstart = i

 #Shower condition true: Steepness Check, Shower Start Layer i or i - 1?

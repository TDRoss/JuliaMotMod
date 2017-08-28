function MotorModelv01(nmpp::Int,nmpw::Int,nmnn::Int,nf::Int,vp::Float64,vw::Float64,vn::Float64,flen::Float64,motln::Float64,boundlen::Float64,koff::Float64,kon::Float64)
#nmpp = number of positive-positive motor pairs
#nmpw = number of positive-wild motor pairs
#nmnn = number of negative-negative motor pairs
#nf = number of filaments
#vp = velocity of positive motors
#vw = velocity of wild motors
#vn = velocity of negative motors
#flen = filament length
#motln = motor length
#boundlen = length of boundary
#koff = motor-filament off rate
#kon = motor-filament on rate

  #Initialize costants
  dt = 10^-4; #time step size (seconds)
  totaltime = 150  #Duration of simulation  (seconds)
  iters = dt:dt:totaltime #range of simulation

  fs = 5 #stall force (pN)
  ks = 0.01 #spring constant (pN/nm)
  fdragll = 10^ -5 #parallel drag on filament
  fdragt = 2*10^-5 #orthogonal drag on filament
  fdragr = 2*10^2 # rotational drag on filament
  hflng = flng/2; #half the filament length
  filtdiff = sqrt(2*dt*10^5) #translational filament diffusion coefficent
  filrdiff = sqrt(2*dt*10&-2) #rotational filament diffusion coefficent
  mtdiff = sqrt(2*dt*4*10^6) #motor translational diffusion coefficent
  offrate = dt*koff
  onrate = dt*kon

  nm = nmpp+nmpw+nmnn #total number of motors in the system
  maxmotloc = boundlen + flen #upper boundary limit on motor position
  minmotloc = -maxmotloc #lower boundary limit on motor position
  timesaveresolution = 1 #time resolution of output array (seconds)
  indxstep = Int(round(timesaveresolution/dt)) #index for when to save history
  totalindx = Int(floor(length(iters)/indxstep))
  filbound = boundlen*2 #Boundary for filament generation
  motbound = (boundlen + flen)*2 #Boundary for motor generation


  #Initialize filaments
  filID = collect(1:nf) #filament ID
  filcxp = rand(nf)*filbound -filbound/2 #x coordinate of filament center
  filcyp = rand(nf)*filbound -filbound/2 #y coordinate of filament center
  filori = rand(nf)*2*pi #Orientation of filament
  filsxp = similar(filcxp) #x coordinate start of filament
  filexp = similar(filcxp) #x coordinate end of filament
  filsyp = similar(filcyp) #y coordinate start of filament
  fileyp = similar(filcyp) #y coordinate end of filament

  @inbounds for j in 1:nf #Calculate filament coordinates
    filsxp[j] = filcxp[j] - hflng * cos(filori[j])
    filsyp[j] = filcyp[j] - hflng * sin(filori[j])
    filexp[j] = filcxp[j] + flng * cos(filori[j])
    fileyp[j] = filcyp[j] + flnd * sin(filori[j])
  end

  #Initalize motors
  mbfID = zeros(Int,nm) #ID of filament that motor is bound to, 0 if none
  mpairID = zeros(Int,nm) #ID of motor pair, 0 if not in pair
  mID = collect(1:nm) #motor ID
  mfree = ones(Int,nm) #motors that are not bound to a filament
  masingle = zeros(Int,nm) #motors that are active for dimerizing
  numpairs = Int(nm/2) #number of possible motor pairs
  pairlistlogical = ones(Int,numpairs) #Logical index of available motor pair IDs
  pairlistdex = collect(1:numpairs) #motor pair IDs

  mpairs = zeros(Int,npairs,4) #Index of motor pairs, index 3 and 4 point to motor mbfID
  availablepairID = 1 #stores the next available pair ID value

  mxp = rand(nm)*motbound - motbound/2 #x coordinate of motor
  myp = rand(nm)*motbound - motbound/2 #y coordinate of motor
  mdfc = similar(mxp) #distance from bound filament center
  moptype = similar(mxp) #opto type of motor, determines which motor pairs form
  mspeed = similar(mxp) #speed of motor, sign indicates plus or minus directed
    @inbounds for j in 1:nm
      if j<=nmpp
        moptype[j] = 1+iseven(j)
        mspeed[j] = vp*dt
      elseif j<=(nmpp + nmpw)
        moptype[j] = 3+iseven(j)
        mspeed[j] = dt*vp*iseven(j)+vw*isodd(j)*dt
      elseif j<=nm
        moptype[j] = 5+iseven(j)
        mspeed[j] = vn*dt
      end
    end
  moptype[1:nmpp] = repeat([1,2],inner=Int(nmpp/2))
  moptype[(1+nmpp):(nmpp+nmpw)] = repeat([3,4],inner=Int(nmpw/2))
  moptype[(1+nmpp):(nmpp+nmpw)] = repeat([3,4],inner=Int(nmpw/2))


  #Store histories of motor and filament positions
  filhist = zeros(totalindx,nf,4)
  mothist = zeros(totalindx,nm,5)

  #For calculating distance between motor pairs
  paironeID = 0
  pairtwoID = 0
  sx = 0.0
  sy = 0.0
  mpdist = 0.0
  atv = 0.0
  sa1 = 0.0
  sa2 = 0.0


  pairdiff = 0.0 #For calculating diffusion on free motor pairs
  boundone = 0 #holds the id of the bound motor in a motor pair
  mpos = zeros(Float64,2) #motor position for querrying KDTree
  sfx = zeros(Float64,nf,2) #sorted filament x start and end point
  sfy = zeros(Float64,nf,2) #sorted filament y start and end point
  #For calulating force on filaments
  forceonfils = similar(filcxp)
  #Begin simulation
  @inbounds for t in iters
  #Calculate net forces on filaments due to motor pair springs
  forceonfils = zeros(nf)
  #Motor pairs that are filament bound
  #requires that 1) motors are paired and 2) each motor is bound
  #we may be able to make a list of both seperately as a set of boolean arrays
  #where each condition is met
    # @inbounds for j in 1:npairs
    #   if mpairs[j,3] >0 && mpairs[j,4] > 0
    #     paironeID = mpairs[j,1]
    #     pairtwoID = mpairs[j,2]
    #     sx = mxp[paironeID] - mxp[pairtwoID]
    #     sy = myp[paironeID] - myp[pairtwoID]
    #     l = sqrt(sx^2+sy^2)
    #     atv = atan2(sy,sx)
    #     sa1 = atv-th1
    #     sa2 = atv-th2 + pi
    #     forceonfils[mpairs[]]
    #   end
    # end
  #Apply all forces to filaments

  #Break up filaments into points for searching
  @inbounds for j in 1:nf
    sfx[j,:] = sort([filsxp[j],filexp[j]])
    sfy[j,:] = sort([filsyp[j],fileyp[j]])
  end

  #Apply forces on motors
  @inbounds for j in 1:nm
    if mpairID[j] == 0 #single motors
      if mbfID[j] == 0 #not bound to filament
        mxp[j] += randn()*mtdiff
        myp[j] += randn()*mtdiff
        #reflect motor back if out of bounds
        if abs(mxp[j]) > maxmotloc
          mxp[j] = maxmotloc*2.0 - mxp[j]
        end
        if abs(myp[j]) > maxmotloc
          myp[j] = maxmotloc*2.0 - myp[j]
        end
      else #bound to filament
        mdfc[j]+= mspeed[j]
        mxp[j] = mdfc[j]*cos(filori[mbfID[j]])+filcxp[mbfID[j]]
        myp[j] = mdfc[j]*sin(filori[mbfID[j]])+filcyp[mbfID[j]]
      end
    else #paired motor
      if mpairs[mpairID[j],3] == 0 && mpairs[mpairID[j],4] == 0 && mpairs[mpairID[j],1] == j #both motors not bound to filament
        pairtwoID = view(mpairs,mpairID[j],2)
        pairdiff = randn()*mtdiff
        mxp[j]+= pairdiff
        mxp[pairtwoID]+= pairdiff
        pairdiff = randn()*mtdiff
        myp[j]+= pairdiff
        myp[pairtwoID]+= pairdiff
        #reflect out of bound motor pairs
        if abs(mxp[j]) > maxmotloc
          mxp[j] = maxmotloc*2.0 - mxp[j]
        end
        if abs(myp[j]) > maxmotloc
          myp[j] = maxmotloc*2.0 - myp[j]
        end
      elseif xor(mpairs[mpairID[j],3] == 0, mpairs[mpairID[j],4] == 0) && mpairs[mpairID[j],1] == j #one motor is bound
        paironeID = view(mpairs,mpairID[j],1)
        pairtwoID = view(mpairs,mpairID[j],2)
        boundone = mpairs[mpairID[j],1] * (mpairs[mpairID[j],3]) + mpairs[mpairID[j],2] * (mpairs[mpairID[j],4])
        mdfc[boundone] += mspeed[boundone]
        pairdiff = mdfc[boundone]*cos(filori[mbfID[boundone]])+filcxp[mbfID[boundone]]-mxp[boundone]
        mxp[paironeID]+= pairdiff
        mxp[pairtwoID]+= pairdiff
        pairdiff = mdfc[boundone]*sin(filori[mbfID[boundone]])+filcyp[mbfID[boundone]]-myp[boundone]
        myp[paironeID]+= pairdiff
        myp[pairtwoID]+= pairdiff
      elseif   pairlistlogical(mpairID[j]) #both motors in pair bound to filament
        mdfc[j]+= mspeed[j]
        mxp[j] = mdfc[j]*cos(filori[mbfID[j]])+filcxp[mbfID[j]]
        myp[j] = mdfc[j]*sin(filori[mbfID[j]])+filcyp[mbfID[j]]
      end
    end

    activetree = BruteTree([mxp[mID[masingle]] myp[mid[masingle]]]) #KDTree of active single motor
    sfx
    @inbounds for j in 1:nm
      if mfree[j] #Find free motors for filament binding
        if mpairID[j] == 0
        checkxp = mxp[j]
        checkyp = myp[j]
        elseif mpairs[mpairID[j],3] == 0 && mpairs[mpairID[j],4] == 0
        paironeID = mpairs[mpairID[j],1]
        pairtwoID = mpairs[mpairID[j],2]
        checkxp = (mxp[paironeID]+mxp[pairtwoID])/2.0
        checkyp = (myp[paironeID]+myp[pairtwoID])/2.0
        else
          if mpairs[mpairID[j],3] == 1
            paironeID = mpairs[mpairID[j],1]
          else
            paironeID = mpairs[mpairID[j],2]
          end
        checkxp = mxp[paironeID]
        checkyp = myp[paironeID]
        end
        @inbounds for p in 1:nf#check if filament is a potential partner
          if sfx[p,1]<checkxp<sfx[p,2] && sfy[p,1]<checkyp[j]<sfy[p,2] #Check if filament is in range
            [t1,t2]=circleintline(checkxp,checkyp,filsxp[p],filsyp[p],flng,filexp[p],fileyp[p],100)
            if 1> t1 > 0  #if filament is in range
              if rand()<onrate
                mdfc[j] = t1*flng
                mxp[j] = filsxp[p] + mdfc[j]*cos(filori)
                myp[j] = filsyp[p] + mdfc[j]*sin(filori)
                mdfc[j] -= hflng
                mfree[j] = 0
                mbfID[j] = p
                break
              end
            elseif 1>t2>0
              if rand()<onrate
                mdfc[j] = t1*flng
                mxp[j] = filsxp[p] + mdfc[j]*cos(filori)
                myp[j] = filsyp[p] + mdfc[j]*sin(filori)
                mdfc[j] -= hflng
                mfree[j] = 0
                mbfID[j] = p
                break
              end
            end
          end
        end
      end
      if !mfree[j] #motors that are filament bound
        if mdfc >= hflng || (rand() < offrate) #motors at end of filament fall off
          mfree[j] = 1
          mbfID[j] = 0
        end
      end

      if masingle[j] == 1 && mpairID[j] == 0  #Find active motors that are available for binding
        mpos = [mxp[j],myp[j]] #holds coordinate of current motor
        idxmx = inrange(activetree,mpos,110,false) #find index of nearby points at max radius
        idxmn = inrange(activetree,mpos,90,false) #find index of nearby points at min radius
        idxs = intersect(idxmx,idxmn) #possible indices
        @inbounds for i in idxs #loop over possible pair partners
          if mpairID[mID[masingle][i]] == 0 #If potential motor is not already paired, pair motors
            masingle[j] = 0 #set motors no longer single
            masingle[mID[masingle][i]] = 0
            #assign pair ID
            mpairID[j] = pairlistdex[pairlistlogical][1]
            mpairID[mID[masingle][i]] == mpairID[j]
            #Remove pair from available ID list
            mairlistlogical[mpairID[j]] = 0
            #Update mpairs vector
            mpairs[mpairID[j],1] = j
            mpairs[mpairID[j],2] = mID[masingle][i]
            mpairs[mpairID[j],3] = view(mbfID,j)
            mpairs[mpairID[j],4] = view(mbfID,mID[masingle][i])
            break
          end
        end
    end
  #Make activated motors dimerize

  #Make unactivated motors monomerize

  #Bind free motors to filaments

  #Detatch bound motors from filaments

  end
  #

  function circleintline(xm::Float64,ym::Float64,xs::Float64,ys::Float64,l::Float64,xe::Float64,ye::Float64,r::Float64)
  a = l^2
  b = 2*(xe-xs)*(xs-xm) + 2*(ye-ys)*(ys-ym)
  c = (xs - xm)^2 + (ys - ym)^2 - r^2

  t1 = 2*c/(-b+sqrt(b^2-4*a*c))
  t2 = 2*c/(-b-sqrt(b^2-4*a*c))
  return t1,t2
  end

using MAT
using NearestNeighbors
function MotorModelv01(outname::String,nmpp::Int,nmpw::Int,nmnn::Int,nf::Int,vp::Float64,vw::Float64,vn::Float64,flen::Float64,motln::Float64,boundlen::Float64,koff::Float64,kon::Float64,preform::Bool)
#outname = string for output file
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
  dt = 0.5*10^(-5.); #time step size (seconds)
  skiptime = Int(floor(5*10^(-2.)/dt)) #time between collision detection events
  totaltime = 20.0  #Duration of simulation  (seconds)
  iters = dt:dt:totaltime #range of simulation

  fs = 5 #stall force (pN)
  ks = 0.1 #spring constant (pN/nm)
  fdragll = 10^-5. #parallel drag on filament
  fdragt = 2*10^-5. #orthogonal drag on filament
  fdragr = 2*10^2. # rotational drag on filament
  hflen = flen/2.0; #half the filament length
  filtdiff = sqrt(2*dt*2*10^5.) #translational filament diffusion coefficent
  filrdiff = sqrt(2*dt*2*10^-2.) #rotational filament diffusion coefficent
  mtdiff = sqrt(2*dt*4*10^6.) #motor translational diffusion coefficent
  offrate = skiptime*dt*koff
  onrate = skiptime*dt*kon

  nm = nmpp+nmpw+nmnn #total number of motors in the system
  maxmotloc = boundlen + flen #upper boundary limit on motor position
  minmotloc = -maxmotloc #lower boundary limit on motor position

  filbound = boundlen*2 #Boundary for filament generation
  motbound = (boundlen + flen)*2#Boundary for motor generation


  #Initialize filaments
  filID = collect(1:nf) #filament ID
  filcxp = rand(nf)*filbound -filbound/2 #x coordinate of filament center
  filcyp = rand(nf)*filbound -filbound/2 #y coordinate of filament center
  filori = rand(nf)*2*pi #Orientation of filament
  filsxp = similar(filcxp) #x coordinate start of filament
  filexp = similar(filcxp) #x coordinate end of filament
  filsyp = similar(filcyp) #y coordinate start of filament
  fileyp = similar(filcyp) #y coordinate end of filament
  filcos = similar(filcxp) #cos of filament angle
  filsin = similar(filcyp) #sin of filament angle

  # For testing two filaments
  # filcxp = [0.0, 0.0]
  # filcyp = [0.0, 0.0]
  # filori = [0.0, pi/2.0]
  @inbounds for j in 1:nf #Calculate filament coordinates
    filcos[j] = cos(filori[j])
    filsin[j] = sin(filori[j])
    filsxp[j] = filcxp[j] - hflen * filcos[j]
    filsyp[j] = filcyp[j] - hflen * filsin[j]
    filexp[j] = filcxp[j] + flen * filcos[j]
    fileyp[j] = filcyp[j] + flen * filsin[j]
  end

  #Initalize motors
  mbfID = zeros(Int,nm) #ID of filament that motor is bound to, 0 if none
  mpairID = zeros(Int,nm) #ID of motor pair, 0 if not in pair
  mID = collect(1:nm) #motor ID
  mfree = ones(Bool,nm) #motors that are not bound to a filament
  masingle = ones(Bool,nm) #motors that are active for dimerizing
  numpairs = Int(nm/2) #number of possible motor pairs
  pairlistlogical = ones(Bool,numpairs) #Logical index of available motor pair IDs
  pairlistdex = collect(1:numpairs) #motor pair IDs

  mpairs = zeros(Int,numpairs,4) #Index of motor pairs, index 3 and 4 point to motor mbfID
  activetimer = zeros(Int,nm) #stores remaining time of motor pairing activation
  timeactive = Int(3) #activation length for motor pairing
  mxp = rand(nm)*motbound - motbound/2 #x coordinate of motor
  myp = rand(nm)*motbound - motbound/2 #y coordinate of motor
  mdfc = similar(mxp) #distance from bound filament center
  moptype = similar(activetimer) #opto type of motor, determines which motor pairs form
  mspeed = similar(mxp) #speed of motor, sign indicates plus or minus directed

  #make motors start out in pairs
  pc = 1
  @inbounds for j in 1:nm #assign parameters to motors
    if preform
      #make motors start out in pairs
      mpairs[pc,1+iseven(j)] = j
      mpairID[j] = pc
      masingle[j] = 0
      pairlistlogical[pc] = 0
      pc+=iseven(j)
      if iseven(j)
        mxp[j] = mxp[j-1]+100
        myp[j] = myp[j-1]
      end
    end

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

  # #For testing motor pair on two filaments
  # mxp = [0.0,70.71]
  # myp = [70.71,0.0]
  # mdfc = [70.71,70.71]
  # mpairs[1,:] = [1,2,1,2]
  # pairlistlogical = zeros(Bool,1)
  # masingle = zeros(Bool,2)
  # mfree = zeros(Bool,2)
  # mbfID = [1,2]
  # mpairID =[1,1]

  #Store histories of motor and filament positions
  timesaveresolution = .1 #time resolution of output array (seconds)
  indxstep = Int(round(timesaveresolution/dt)) #index for when to save history
  totalindx = Int(floor(length(iters)/indxstep))+1
  filhist = zeros(totalindx,nf,3) #sxp,syp,ori
  mothist = zeros(totalindx,nm,3) #xp,yp,filament ID
  forcehist = zeros(totalindx,nf,3) #xf,yf,rf last dimension is motor generated/random forces
  steptrigger = Int(round(timesaveresolution/dt))
  tc = 1
  hc = 1
  filhist[hc,:,1] = filsxp
  filhist[hc,:,2] = filsyp
  filhist[hc,:,3]  = filori
  mothist[hc,:,1] = mxp
  mothist[hc,:,2] = myp
  #For calculating distance between motor pairs
  #Motor IDs
  paironeID = 0
  pairtwoID = 0
  #Filament IDs
  filoneID = 0
  filtwoID = 0
  #Distance between motors in x and y
  sx = 0.0
  sy = 0.0
  mpdist = 0.0
  #angle between two motors
  atv = 0.0
  sa1 = 0.0
  sa2 = 0.0
  #spring force on motors
  spforce = 0.0
  #Calculate when motor will fall off filament
  mdcoff = hflen - dt*skiptime*vp
  #holds valid optype for motor pairing
  validoptype = 0
  #holds stochastic forces applied to filaments
  randforce = zeros(Float64, 3)
  #Active region positions
  # actreg = [-5*flen hflen -5*flen 5*flen; flen 5*flen -5*flen 5*flen]  #xs xe ys ye; region 2
  linkreg = [-5*flen 5*flen -hflen hflen]
  actreg = [-10*flen 10*flen -10*flen 10*flen;-10*flen 10*flen -10*flen 10*flen]

  pairdiff = 0.0 #For calculating diffusion on free motor pairs
  boundone = 0 #holds the id of the bound motor in a motor pair
  mpos = zeros(Float64,2) #motor position for querrying KDTree
  sfx = zeros(Float64,nf,2) #sorted filament x start and end point
  sfy = zeros(Float64,nf,2) #sorted filament y start and end point
  #For calulating force on filaments
  forceonfils = similar(filcxp)
  skc = 0
  forceonfils = zeros(Float64,nf,3)
  forceonmots = zeros(Float64,nm)
@inbounds  for t in iters   #Begin simulation
    skc+=1
    tc+=1
    if tc == steptrigger
      hc+=1
      println(t/totaltime)
    end
    forceonfils .= 0.0
    forceonmots .= 0.0

     @inbounds for j in 1:numpairs #Calculate net forces on filaments due to motor pair springs
      if mpairs[j,3] >0 && mpairs[j,4] > 0 #If motor pairs are both filament bound

        #Get motor and filament IDs
        paironeID = mpairs[j,1]
        pairtwoID = mpairs[j,2]
        filoneID = mpairs[j,3]
        filtwoID = mpairs[j,4]
        #Calculate distance and angle between motor pairs
        sx = mxp[pairtwoID] - mxp[paironeID]
        sy = myp[pairtwoID] - myp[paironeID]
        l = sqrt(sx^2+sy^2)
        if l > motln
          atv = atan2(sy,sx)
          sa1 = atv-filori[filoneID]
          sa2 = atv-filori[filtwoID] + pi
          #Calculate spring force on motors
          spforce= ks * (l-motln)

          forceonmots[paironeID] = abs(spforce)
          forceonmots[pairtwoID] =  forceonmots[paironeID]

          #Calculate force applied to filaments
          forceonfils[filoneID,1] += spforce*cos(sa1)
          forceonfils[filtwoID,1] += spforce*cos(sa2)
          forceonfils[filoneID,2] += spforce*sin(sa1)
          forceonfils[filtwoID,2] += spforce*sin(sa2)
          forceonfils[filoneID,3] += mdfc[paironeID] * spforce * sin(sa1)
          forceonfils[filtwoID,3] += mdfc[pairtwoID] * spforce * sin(sa2)

        end
      end
    end


     @inbounds for j in 1:nf #Apply all forces to filaments
      # randforce = randn(3)
      filcxp[j] += dt*(forceonfils[j,1]*filcos[j]/fdragll-forceonfils[j,2]*filsin[j]/fdragt) #+ randn()*filtdiff
      filcyp[j] += dt*(forceonfils[j,1]*filsin[j]/fdragll+forceonfils[j,2]*filcos[j]/fdragt) #+ randn()*filtdiff
      filori[j] += dt*forceonfils[j,3]/fdragr #+ randn()*filrdiff

      #update filament start and end points
      filcos[j] = cos(filori[j])
      filsin[j] = sin(filori[j])
      filsxp[j] = filcxp[j] - hflen * filcos[j]
      filsyp[j] = filcyp[j] - hflen * filsin[j]
      filexp[j] = filcxp[j] + hflen * filcos[j]
      fileyp[j] = filcyp[j] + hflen * filsin[j]


      if tc == steptrigger #update filament history array
        filhist[hc,j,1] = filsxp[j]
        filhist[hc,j,2] = filsyp[j]
        filhist[hc,j,3]  = filori[j]
        forcehist[hc,j,1] = forceonfils[j,1]
        forcehist[hc,j,2] = forceonfils[j,2]
        forcehist[hc,j,3] = forceonfils[j,3]
      end
      if skc == skiptime
        #Get sorted filament start and end points, used for collision detection with motors
        if filcos[j]>0
          sfx[j,1] = filsxp[j] - 100
          sfx[j,2] = filexp[j] + 100
        else
          sfx[j,1] = filexp[j] - 100
          sfx[j,2] = filsxp[j] + 100
        end
        if filsin[j]>0
          sfy[j,1] = filsyp[j] - 100
          sfy[j,2] = fileyp[j] + 100
        else
          sfy[j,1] = fileyp[j] - 100
          sfy[j,2] = filsyp[j] + 100
        end
      end
    end

    #Apply forces on motors
     @inbounds for j in 1:nm
      if mpairID[j] == 0 #single motors
        if mbfID[j] == 0 #not bound to filament
          mxp[j] += randn()*mtdiff
          myp[j] += randn()*mtdiff
          #reflect motor back if out of bounds
          if abs(mxp[j]) > maxmotloc
            mxp[j] = sign(mxp[j])*maxmotloc*2.0 - mxp[j]
          end
          if abs(myp[j]) > maxmotloc
            myp[j] = sign(myp[j])*maxmotloc*2.0 - myp[j]
          end
        else #bound to filament
          mdfc[j]+= mspeed[j]
          mxp[j] = mdfc[j]*filcos[mbfID[j]]+filcxp[mbfID[j]]
          myp[j] = mdfc[j]*filsin[mbfID[j]]+filcyp[mbfID[j]]
        end
      else #paired motor
        if mpairs[mpairID[j],1] == j && mpairs[mpairID[j],3] == 0 && mpairs[mpairID[j],4] == 0 #both motors not bound to filament
          pairtwoID = mpairs[mpairID[j],2]
          pairdiff = randn()*mtdiff
          mxp[j]+= pairdiff
          mxp[pairtwoID]+= pairdiff
          pairdiff = randn()*mtdiff
          myp[j]+= pairdiff
          myp[pairtwoID]+= pairdiff
          #reflect out of bound motor pairs
          if abs(mxp[j]) > maxmotloc
            mpdist = 2.0*(sign(mxp[j])*maxmotloc - mxp[j])
            mxp[j] += mpdist
            mxp[pairtwoID] += mpdist
          end
          if abs(myp[j]) > maxmotloc
            mpdist = 2.0*(sign(myp[j])*maxmotloc - myp[j])
            myp[j] += mpdist
            myp[pairtwoID] += mpdist
          end
        elseif mpairs[mpairID[j],1] == j && xor(mpairs[mpairID[j],3] == 0, mpairs[mpairID[j],4] == 0) #one motor is bound
          paironeID = mpairs[mpairID[j],1]
          pairtwoID = mpairs[mpairID[j],2]
          if mpairs[mpairID[j],3]>0
            boundone = mpairs[mpairID[j],1]
          else
            boundone = mpairs[mpairID[j],2]
          end
          mdfc[boundone] += mspeed[boundone]
          pairdiff = mdfc[boundone]*filcos[mbfID[boundone]]+filcxp[mbfID[boundone]]-mxp[boundone]
          mxp[paironeID]+= pairdiff
          mxp[pairtwoID]+= pairdiff
          pairdiff = mdfc[boundone]*filsin[mbfID[boundone]]+filcyp[mbfID[boundone]]-myp[boundone]
          myp[paironeID]+= pairdiff
          myp[pairtwoID]+= pairdiff
        elseif   mpairs[mpairID[j],3] >0 && mpairs[mpairID[j],4] >0  #both motors in pair bound to filament
          if forceonmots[j] < fs
            mdfc[j]+= mspeed[j]*(1-forceonmots[j]/fs)
          end
          mxp[j] = mdfc[j]*filcos[mbfID[j]]+filcxp[mbfID[j]]
          myp[j] = mdfc[j]*filsin[mbfID[j]]+filcyp[mbfID[j]]
        end
      end
      if tc == steptrigger #update motor history array
        mothist[hc,j,1] = mxp[j]
        mothist[hc,j,2] = myp[j]
        mothist[hc,j,3] = mbfID[j]
      end
    end

      if skc==skiptime
        @inbounds for j in 1:nm
          if activetimer[j] >0 #Count down active timer for motor pairing
            activetimer[j] -= 1
          end
          #Check if motor is in active region
          if (actreg[1,1]<mxp[j]<actreg[1,2] && actreg[1,3]<myp[j]<actreg[1,4]) || (actreg[2,1]<mxp[j]<actreg[2,2] && actreg[2,3]<myp[j]<actreg[2,4]) || (t>140.0 && linkreg[1]<mxp[j]<linkreg[2] && linkreg[3]<myp[j]<linkreg[4])
            activetimer[j] = timeactive
            if mpairID[j] == 0
              masingle[j] = 1
            end
          end
          if activetimer[j] == 0
            if mpairID[j]>0 #detatch deactivated motor pairs
              paironeID = mpairs[mpairID[j],1]
              pairtwoID = mpairs[mpairID[j],2]
              if paironeID == j
                masingle[pairtwoID] = 1
              else
                masingle[paironeID] = 1
              end
              mpairs[mpairID[j],:] = 0
              pairlistlogical[mpairID[j]] = 1
              mpairID[paironeID] = 0
              mpairID[pairtwoID] = 0
            end
            masingle[j] = 0
          end
          if mfree[j] #Find free motors for filament binding
            boundone = 0
            if mpairID[j] == 0 #single motors
              checkxp = mxp[j]
              checkyp = myp[j]
            elseif mpairs[mpairID[j],3] == 0 && mpairs[mpairID[j],4] == 0 #free pair
              paironeID = mpairs[mpairID[j],1]
              pairtwoID = mpairs[mpairID[j],2]
              checkxp = (mxp[paironeID]+mxp[pairtwoID])/2.0
              checkyp = (myp[paironeID]+myp[pairtwoID])/2.0
            else #motor is part of pair with one filament bound
              if mpairs[mpairID[j],3] >0 #Get ID of bound motor, use that as center of binding radius
                paironeID = mpairs[mpairID[j],1]
                boundone = mpairs[mpairID[j],3]
              else
                paironeID = mpairs[mpairID[j],2]
                boundone = mpairs[mpairID[j],4]
              end
            checkxp = mxp[paironeID]
            checkyp = myp[paironeID]
            end
             @inbounds for p in 1:nf#check if filament is a potential partner
              if sfx[p,1]<checkxp<sfx[p,2] && sfy[p,1]<checkyp<sfy[p,2] && boundone != p #Check if filament is in range and exclude motor pairs binding to same filament
                t1,t2=circleintline(checkxp,checkyp,filsxp[p],filsyp[p],flen,filexp[p],fileyp[p],100.0)
                if 1.0> t1 > 0.0  #if filament is in range
                  if rand()<onrate
                    mdfc[j] = t1*flen-hflen
                    mxp[j] = filcxp[p] + mdfc[j]*filcos[p]
                    myp[j] = filcyp[p] + mdfc[j]*filsin[p]
                    mfree[j] = 0
                    mbfID[j] = p
                    if mpairID[j] >0
                      if mpairs[mpairID[j],1] == j
                        filoneID = 3
                      else
                        filoneID = 4
                      end
                      mpairs[mpairID[j],filoneID] = p
                    end
                    break
                  end
                elseif 1.0>t2>0.0
                  if rand()<onrate
                    mdfc[j] = t1*flen-hflen
                    mxp[j] = filcxp[p] + mdfc[j]*filcos[p]
                    myp[j] = filcyp[p] + mdfc[j]*filsin[p]
                    mfree[j] = 0
                    mbfID[j] = p
                    if mpairID[j] >0
                      if mpairs[mpairID[j],1] == j
                        filoneID = 3
                      else
                        filoneID = 4
                      end
                      mpairs[mpairID[j],filoneID] = p
                    end
                    break
                  end
                end
              end
            end
          end

          if !mfree[j] #motors that are filament bound, see if they unbind
            if mdfc[j] >= mdcoff || (rand() < offrate) #motors at end of filament fall off
              mfree[j] = 1
              mbfID[j] = 0
              if mpairID[j] >0
                if mpairs[mpairID[j],1] == j
                  filoneID = 3
                else
                  filoneID = 4
                end
                mpairs[mpairID[j],filoneID] = 0
              end
            elseif rand() < offrate*exp(forceonmots[j]/2.5) #koff of motors from filaments
              mfree[j] = 1
              mbfID[j] = 0
              if mpairID[j] >0
                if mpairs[mpairID[j],1] == j
                  filoneID = 3
                else
                  filoneID = 4
                end
                mpairs[mpairID[j],filoneID] = 0
              end
            end
          end
        end

      if sum(masingle)>0 #make tree if single motors are available for pairing
        activetree = KDTree([mxp[mID[masingle]] myp[mID[masingle]]]') #KDTree of active single motor
        treedex = mID[masingle]
      end

      @inbounds for j in 1:nm #loop over motors for pairing
        if masingle[j] == 1 && mpairID[j] == 0  #Find active motors that are available for pairing
          mpos[1] = mxp[j]
          mpos[2] = myp[j] #holds coordinate of current motor
          if iseven(moptype[j])
            validoptype = moptype[j]-1 #valid optype binding partner
          else
            validoptype = moptype[j]+1
          end
          idxmx = inrange(activetree,mpos,110,false) #find index of nearby points at max radius
          idxmn = inrange(activetree,mpos,90,false) #find index of nearby points at min radius
          idxs = setdiff(idxmx,idxmn) #possible indices
          @inbounds for i in idxs #loop over possible pair partners
            if mpairID[treedex[i]] == 0 && validoptype==moptype[treedex[i]] &&(mbfID[j]==0 || mbfID[treedex[i]]==0 || mbfID[j] != mbfID[treedex[i]]) #If motor is not already paired or on same filament, pair motors
              boundone = treedex[i]
              masingle[j] = 0 #set motors no longer single
              masingle[boundone] = 0
              #assign pair ID
              mpairID[j] = pairlistdex[pairlistlogical][1]

              mpairID[boundone] = mpairID[j]
              #Remove pair from available ID list
              pairlistlogical[mpairID[j]] = 0
              #Update mpairs vector
              mpairs[mpairID[j],1] = j
              mpairs[mpairID[j],2] = boundone
              mpairs[mpairID[j],3] = mbfID[j]
              mpairs[mpairID[j],4] = mbfID[boundone]
              break
            end
          end
        end
      end
    skc = 0
  end

    if tc == steptrigger
      tc = 0
    end
  end
    file = matopen(outname,"w")
    write(file,"filhist",filhist)
    write(file,"mothist",mothist)
    write(file,"forcehist",forcehist)
    close(file)
    # return mothist, filhist
end


  function circleintline(xm::Float64,ym::Float64,xs::Float64,ys::Float64,l::Float64,xe::Float64,ye::Float64,r::Float64)
  a = l^2.0
  b = 2.0*(xe-xs)*(xs-xm) + 2.0*(ye-ys)*(ys-ym)
  c = (xs - xm)^2.0 + (ys - ym)^2.0 - r^2.0
  disc = (b)^2.0-4.0*a*c
  if disc >=0
    t1 = 2.0*c/(-b+sqrt(disc))
    t2 = 2.0*c/(-b-sqrt(disc))
  else
    t1 = -1.0
    t2 = -1.0
  end

  return t1,t2
  end

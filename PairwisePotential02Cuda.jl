using JLD
using MAT
using Distances
using CuArrays

#Pairwise potential function to describe filament interactions
#Difference between version 01 and 02 is that 02 does not have a limited
#interaction radius and instead localizes particle interactions to occur
#within defined boundaries
function pairwisepotential(argvec)
    CuArrays.allowscalar(false)
    outname=argvec[1]::String #Name of output file
    nf = argvec[2]::Int #Number of filaments
    ffdim = argvec[3]::Float32 #Filament field dimension
    flen = argvec[4]::Float32 #Filament length (sets interaction cutoff)
    mspeed = argvec[5]::Float32 #Max Speed of attraction
    # magp  = argvec[6]::Float32 #Magnitude of potential

    #Initialize constants
    #######
    #set  time parameters
    dt = 2.5f0 #time step size
    totaltime = 2400.0f0 #total simulation time
    iters = dt:dt:totaltime
    dspeed = mspeed*dt

    #Calculate diffusion parameters
    avdrag = ((2.0f0*pi*8.9f0*10^(-10.0f0)*flen/(log(flen/(24.0f0))-0.2f0)) +
     (4.0f0*pi*8.9f0*10^(-10.0f0)*flen/(log(flen/(24.0f0))-0.84f0)) )/2 #Average drag on filament
    avdiff = sqrt(2.0f0*dt*1.0f0/avdrag) #Diffusion on filament

    #Initialize filaments
    # filpos = cu(zeros(Float32,2,nf)) #Filament positions x,y
    # filpos = cu(rand(Float32,2,nf))*ffdim*2 -ffdim #Filament  position
    filx = cu(rand(Float32,nf))*ffdim*2 -ffdim
    fily = cu(rand(Float32,nf))*ffdim*2 -ffdim



    #History of positions
    timesaveresolution = 10.0f0 #time resolution of output array (seconds)
    indxstep = Int(round(timesaveresolution/dt)) #index for when to save history
    totalindx = Int(floor(length(iters)/indxstep))+1
    filhist = (zeros(Float32,totalindx,2,nf)) #xp, yp

    xx = ([i for i in -40000:500:40000, j in -40000:500:40000][:])
    yy = ([j for i in -40000:500:40000, j in -40000:500:40000][:])
    mgrid = transpose([xx yy])
    pothist = zeros(totalindx,length(xx))

    #Save first time step
    hc = 1 #Index of saved states
    tc = 1 #Time past between saved states
    filhist[hc,1,:] .= collect(filx)
    filhist[hc,2,:] .= collect(fily)
    # dif = zeros(2)
    magc = 1.0f0*10.0f0^4.0f0

    @inbounds for t in iters   #Begin simulation
        tc+=1

        #Calculate diffusion forces on filaments
        filfx = cu(randn(nf)).*avdiff
        filfy = cu(randn(nf)).*avdiff
        # filforce = zeros(2,nf)

        #Identify filaments in excitation zone
        vindx = isinzone(filx,fily,t)

        #Calculate pairwise distance
        # pwdist = cu(pairwise(Euclidean(),collect(filpos[:,indxisinzone])))

        #get distances in interaction range
        # grng = 0.<pwdist
        # ninter = sum(grng[i,:])#number of interactions
        # diff = cu(zeros(2,length(indxisinzone),length(indxisinzone)))

        diffx = broadcast(-,transpose(filx),filx)
        diffy = broadcast(-,transpose(fily),fily)
        # pwdist = sqrt.(diffx.^2 .+ diffy.^2)
        # inzn = vindx.*(0.<pwdist.<flen).*transpose(vindx)

        # normy = sum(inzn,2)[:,1]
        # normy += (normy.<1)
        theta = CUDAnative.atan2.(diffy,diffx)
        ct = CUDAnative.cos.(theta)
        st = CUDAnative.sin.(theta)
        # mag = dspeed./(1.0f0.+exp.(pwdist.-flen))
        for i = 1:3
        inzn = broadcast(&,vindx .== i,transpose(vindx .== i))
        normy = sum(inzn,2)[:,1]
        normy += (normy.==0)
        vx = dspeed.*sum(inzn.*ct,2)[:,1]./normy
        vy = dspeed.*sum(inzn.*st,2)[:,1]./normy
        filfx .+= vx
        filfy .+= vy
        end


        # #Calculate pairwise attraction
        # @inbounds for i in 1:length(indxisinzone)
        #     ninter = sum(grng[i,:])#number of interactions
        #     if ninter>0
        #         diff = filpos[:,indxisinzone[grng[i,:]]] .- filpos[:,indxisinzone[i]]
        #         theta = CUDAnative.atan2.(diff[2,:],diff[1,:])
        #         mag = dspeed./(1.+exp.(pwdist[i,find(grng[i,:])]-flen))
        #
        #         vx = mean(mag.*CUDAnative.cos.(theta))
        #         vy = mean(mag.*CUDAnative.sin.(theta))
        #
        #         filforce[1,indxisinzone[i]]+= vx
        #         filforce[2,indxisinzone[i]]+= vy
        #
        #     end
        # end

        filx .+= filfx
        fily .+= filfy

        #reflect filament back if out of bounds
        w = abs.(filx) .> ffdim
        nw = .~w
        filx .= w.*(sign.(filx).*ffdim.*2.0f0 .- filx) .+ nw.*filx
        w .= abs.(fily) .> ffdim
        nw .= .~w
        fily .= w.*(sign.(fily).*ffdim.*2.0f0 .- fily) .+ nw.*fily

        # @inbounds for j  in 1:nf
        #     if abs(filpos[1,j]) > ffdim
        #       filpos[1,j] = sign(filpos[1,j])*ffdim*2.0f0 - filpos[1,j]
        #     end
        #     if abs(filpos[2,j]) > ffdim
        #       filpos[2,j] = sign(filpos[2,j])*ffdim*2.0f0 - filpos[2,j]
        #     end
        # end

        if tc == indxstep #check if it is time to save state
          hc+=1
          tc=0
          filhist[hc,1,:] .= collect(filx)
          filhist[hc,2,:] .= collect(fily)
          # if hc==2 || hc == round(252/timesaveresolution) || hc == round(1402/timesaveresolution)
          ptinzn = (isinzone(xx,yy,t))
          # end
          for j=1:3
              rinzn = collect(vindx.==j)
              vinzn = ptinzn.==j
              potdist = pairwise(Euclidean(),filhist[hc,:,rinzn],mgrid[:,vinzn])
              # potdist[potdist.>flen] = 0
              pothist[hc,vinzn] = dspeed.*mean(potdist,1)
          end
          # potdist = filx
          # potdist = pairwise(Euclidean(),filhist[hc,:,collect(indxisinzone)],mgrid[:,ptinzn])
          # potdist[potdist.>30000]=0
          # pothist[hc,ptinzn] = mean(potdist,1)
        end


    end
    file = matopen(outname,"w")
    write(file,"filhist",filhist)
    write(file,"argvec",argvec)
    write(file,"pothist",pothist)
    write(file,"xx",xx)
    write(file,"yy",yy)
    close(file)

    # save(string(outname,".jld"), "filhist", filhist, "argvec", argvec)
end

function isinzone(filx,fily,t)
    # indx = ((-30000.0f0.<filx.<-15000.0f0) .& (-30000.0f0.<fily.<30000.0f0)) .|
    # ((-30000.0f0.<filx.<30000.0f0) .& (15000.0f0.<fily.<30000.0f0)) .|
    # ((15000.0f0.<filx.<30000.0f0).&(-30000.0f0.<fily.<30000.0f0)) .|
    # ((-30000.0f0.<filx.<30000.0f0).&(-30000.0f0.<fily.<-15000.0f0))
    # indx = ((-10000.0f0.<filx.<30000.0f0) .& (-10000.0f0.<fily.<10000.0f0)) .|
    #     ((-30000.0f0.<filx.<-10000.0f0) .& (-30000.0f0.<fily.<30000.0f0))

    # if t<250
    # indx = (sqrt.((filx-25000.0f0).^2+(fily+25000.0f0).^2) .<= 10000.0f0) .|
    # (sqrt.((filx).^2+(fily-17000.0f0).^2) .<= 10000.0f0)
    # elseif t<1400
    #     indx =
    #         (sqrt.((filx).^2+(fily-17000.0f0).^2) .<= 10000.0f0) .|
    #     ((-27000.0f0.<filx.<27000.0f0) .& (-27500.0f0.<fily.<-22500.0f0))
    # else
    #     # indx = (sqrt.((filx).^2+(fily-15000.0f0).^2) .<= 10000.0f0) .|
    #     indx = (-30000.0f0.<fily.<19000.0f0) .& (-2500.0f0.<filx.<2500.0f0)
    # end


        if t<250
        indx = 1.*(sqrt.((filx-25000.0f0).^2+(fily+25000.0f0).^2) .<= 10000.0f0) .| 2.*(sqrt.((filx+25000.0f0).^2+(fily+25000.0f0).^2) .<= 10000.0f0) .|
        3.*(sqrt.((filx).^2+(fily-17000.0f0).^2) .<= 10000.0f0)
    elseif t<1400
            indx =    1.*(sqrt.((filx).^2+(fily-17000.0f0).^2) .<= 10000.0f0) .|
            2.*((-27000.0f0.<filx.<27000.0f0) .& (-27500.0f0.<fily.<-22500.0f0))
        else
            # indx = (sqrt.((filx).^2+(fily-15000.0f0).^2) .<= 10000.0f0) .|
            indx = 1.*((-30000.0f0.<fily.<19000.0f0) .& (-2500.0f0.<filx.<2500.0f0))
        end

    # indx = 15000.0f0.<=sqrt.(sum(filpos.^2,1)).<= 30000.0f0
    return indx
end

#Continuous callback

#############################
#Callback 1: healthy → Infect
#############################

#Cells are infected after a certain viral load is met (with μ=800 virons, σ=200 virons)
const maxViralLoad = rand(Normal(m2c(800.0),m2c(200.0)),N,N)

function CheckInfect(u,t,integrator)
    #Reshape the array into 3D
    uReshaped =  reshape(u,N,N,species)

    #Check the whole virus grid for cells with viral conc. greater than 0.4 nM
    return any( @. (uReshaped[:,:,end] > maxViralLoad) & integrator.p.infectFirstAttempt)
end


#These don't need to be constructed every loop
    const Index = CartesianIndices(zeros(N,N))
    #Get the first and last index
    const Ifirst, Ilast = first(Index), last(Index)

function TryInfect(integrator)
    #Probability that an infected cell will spread infection
    chanceOfInfection = Bernoulli(0.5)
    cellsInfected = integrator.p.cellsInfected
    @show integrator.t

    #Sundials does not presereve problem shape...
    uReshaped = reshape(integrator.u,N,N,species)


    #Loop through all of the infectious cells
     for I in findall(uReshaped[:,:,end] .> maxViralLoad)
        #Loop through all neighbors to see if they get infected
        for J in max(Ifirst, I-Ifirst):min(Ilast, I+Ifirst)
            #If the cell is healthy and not the current cell try to infect
            @inbounds if isinf(cellsInfected[J]) && (I != J) && rand(chanceOfInfection)
                #Add Viral DNA to the cell
                uReshaped[J,2] = rand(Uniform(0.0,m2c(1e3)))
                #Mark that the cell was infected
                cellsInfected[J] = integrator.t
            end
        end
        #This cell has used its one attempt to infect neighboring cells
        integrator.p.infectFirstAttempt[I] = false
    end

    integrator.u = uReshaped[:]
end

cbInfect = DiscreteCallback(CheckInfect,TryInfect, save_positions=(false,false))


#############################
#Callback 2: Infect → Dead
#############################

#Cells die some time after infection (with μ=8.0 hours, σ=1.0 hours)
const time2Death = rand(Normal(8.0,1.0),N,N)

function DeathCheck(u,t,integrator)
    cellsInfected = integrator.p.cellsInfected
    #Loop through and check if a cell should die
    cellCounter = 1
    while cellCounter < length(cellsInfected)
        #If cell infected and cell past time of death
        pastTimeOfDeath = (t - cellsInfected[cellCounter]) > time2Death[cellCounter]
        isAlive = isinf(integrator.p.cellsDead[cellCounter])
        if pastTimeOfDeath & isAlive
            return true
        end
        cellCounter += 1
    end
    return false
end

function KillCell(integrator)
    cellsInfected = integrator.p.cellsInfected
    #Find all the cells to kill
    pastTimeOfDeath = (integrator.t .- cellsInfected) .> time2Death
    isAlive = isinf.(integrator.p.cellsDead)
    targets = findall(isAlive .& pastTimeOfDeath)
    #println(sum(integrator.p.deathParameter .== 0.0))
    #Set their deathParameter to stop all cell function
    integrator.p.deathParameter[targets] .= 0.0
    #Mark time of death
    integrator.p.cellsDead[targets] .= integrator.t
end

cbDead = DiscreteCallback(DeathCheck,KillCell,save_positions=(false,false))

cb =CallbackSet(cbInfect,cbDead)

# some annotations as reminders for myself

## Adding MatFileHandler

in the terminal where the `.sln` file is:

    $  dotnet add KTzV2.csproj package MatFileHandler

## Test command lines

 * Running Poisson stimulus test for fixed `J` and variable `r` (Poisson rate) -- stimulus-response function
```
    $ ./KTzV2.exe -run netType=SquareLatticeFreeBC Lx=20 Ly=20 netDir=No dim=2 K=0.6 T=0.35 d=0.001 l=0.008 xR=-0.7 H=0.0 neuron=KTz x0=-0.697156411891772 y0=-0.697156411891772 z0=-0.022748704865822 iCond=ProgramSpecified sType=KTChemicalSynapse tauf=2 taug=2 R=0.036 noiseRatio=0.1 coupParam=Homogeneous stimType=PoissonProcess indChoice=MultipleRandom sampFrac=1.0 nStart=0 nSteps=10000 cVar=NumberOfSpikes samp=Full simType=Bifurcation oFileFormat=mat outAvgMode=OverTime wData=No wCSV=No wAvalDist=No wObs=No nSim=1 bifWrite=InTheEnd writeRhoTS=No dynType=ContinuousTime rRange="logspace(0.000001:1.0:10)" JRange="linspace(-0.12:-0.12:1)" saveSpikeTimes=Yes simTimeScheme=ProportionalToPoissonRate
```

 * Running standard density of active sites phase transition vs. `J`
```
    $ ./KTzV2.exe -run netType=SquareLatticeFreeBC Lx=20 Ly=20 netDir=No dim=2 K=0.6 T=0.35 d=0.001 l=0.008 xR=-0.7 H=0.0 neuron=KTz x0=-0.697156411891772 y0=-0.697156411891772 z0=-0.022748704865822 iCond=ProgramSpecified sType=KTChemicalSynapse tauf=2 taug=2 R=0.036 noiseRatio=0.1 coupParam=Homogeneous stimType=DeltaWhenInactive tBin=50 indChoice=SquareCenter sampFrac=1.0 nStart=0 nSteps=10000 cVar=NumberOfSpikes samp=Full simType=Bifurcation oFileFormat=mat outAvgMode=OverTime wData=No wCSV=No wAvalDist=No wObs=No nSim=1 bifWrite=InTheEnd writeRhoTS=Yes dynType=ContinuousTime IRange="linspace(0.3:0.3:1)" JRange="linspace(0.005:0.020:10)" saveSpikeTimes=Yes
```

## add keyboard shortcut to run task for build

You can define a keyboard shortcut for any task. From the Command Palette (*Ctrl+Shift+P*), select *Preferences: Open Keyboard Shortcuts (JSON)*, 

    {
        "key": "ctrl+shift+b",
        "command": "workbench.action.tasks.runTask",
    }

## install .NET on Ubuntu

https://learn.microsoft.com/en-us/dotnet/core/install/linux-ubuntu#register-the-microsoft-package-repository

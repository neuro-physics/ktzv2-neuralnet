# some annotations as reminders for myself

## Adding MatFileHandler

in the terminal where the `.sln` file is:

    $  dotnet add KTzV2.csproj package MatFileHandler

## Test command lines

```
    $ KTzV2.exe -run netType=SquareLatticeFreeBC Lx=4 Ly=4 netDir=No dim=2 K=0.6 T=0.35 d=0.001 l=0.008 xR=-0.7 H=0 neuron=KTz x0=-0.697156411891772 y0=-0.697156411891772 z0=-0.0227487048658225 iCond=ProgramSpecified sType=KTNoisyChemicalSynapse tauf=2 taug=2 R=0.036 noiseRatio=0.1 coupParam=Homogeneous stimType=PoissonProcess indChoice=Random sampFrac=1 nStart=0 nSteps=10000 cVar=NumberOfSpikes samp=Full simType=Dynamics oFileFormat=mat outAvgMode=OverTime wData=Yes wCSV=No wSpk=No wObs=No nSim=1 bifWrite=OnTheEnd writeRhoTS=Yes dynType=ContinuousTime
```

optional `JRange=linspace(-0.15:-0.25:4)`;    `rRange=logspace(0.01:0.1:4)`;    `dDisorder=Gaussian:0.001:0.0001 TDisorder=Uniform:0.2:0.5`

```
    KTzV2.exe -run netType=SquareLatticeFreeBC Lx=4 Ly=4 netDir=No dim=2 K=0.6 T=0.35 d=0.001 l=0.008 xR=-0.7 H=0 neuron=KTz x0=-0.697156411891772 y0=-0.697156411891772 z0=-0.0227487048658225 iCond=ProgramSpecified sType=KTChemicalSynapse tauf=2 taug=2 R=0.036 noiseRatio=0.1 coupParam=Homogeneous stimType=PoissonProcess indChoice=Random sampFrac=1 nStart=0 nSteps=10000 cVar=NumberOfSpikes samp=Full simType=Dynamics oFileFormat=mat outAvgMode=OverTime wData=Yes wCSV=No wSpk=No wObs=No nSim=1 bifWrite=OnTheEnd writeRhoTS=Yes dynType=ContinuousTime J=-0.2 nJ=1
```

## add keyboard shortcut to run task for build

You can define a keyboard shortcut for any task. From the Command Palette (*Ctrl+Shift+P*), select *Preferences: Open Keyboard Shortcuts (JSON)*, 

    {
        "key": "ctrl+shift+b",
        "command": "workbench.action.tasks.runTask",
    }

## install .NET on Ubuntu

https://learn.microsoft.com/en-us/dotnet/core/install/linux-ubuntu#register-the-microsoft-package-repository

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
//using System.Threading.Tasks;
using KTzV2.Neurons;
using KTzV2.Synapses;
using KTzV2.Stimuli;
using KTzV2.Data.Header;
using KTzV2.Maths.Random;
using KTzV2.Maths.Matrices;
using KTzV2.Maths.Matrices.AdjacencyMatrix;
using KTzV2.Data;
using MatFileHandler;

namespace KTzV2.Sims.Network
{
    public enum StimulusIndexChoice
    {
        Fixed,
        SquareCenter,
        Random,
        MultipleRandom,
        MostConnected
    }

    public enum SimulationTimeScheme
    {
        // total simulation time (nSteps - nStartStep; discounting transient)
        // is proportional to 10/(r*N)
        // r->Poisson rate, N->number of neurons in the network
        ProportionalToPoissonRate,
        // total simulation time is whatever is set by the user
        Free
    }

    public enum CouplingParam
    {
        Homogeneous, // same J for every synapse
        AdjMatrix // uses adj matrix weight as each J
    }

    public class KTzNetworkSimulator
    {
        private SparseMatrix<Double> AMatrix { get; set; }
        private Int32 Lx, Ly, Lz, nSteps, nStartStep;
        private InitialConditionType icType;
        private AdjacencyMatrixType networkType;
        private SynapseType synapseType;
        private StimulusType stimulusType;
        private NeuronType neuronType;
        private NoiseType noiseType;
        private CouplingParam couplingParamType;
        private Boolean doWriteDif = false;
        private Int32 nNeighbours;
        public Int32 nNeurons { get; private set; }
        private Int32 netDim;
        private Int32 tForStimulus;
        private Double x0i, y0i, z0i;
        private Double[] x0p, y0p, z0p;
        private Double K, T, d, l, xR, H, J, tauf, taug, Theta, Q;
        private Double noiseAmp, noiseRatio;
        //private Double[] Ji;
        private GetICInstance icGetter;
        //private Int32[] tData;
        private Double[][] xData;
        private List<Int32> spikeNeuronData;
        private List<Int32> spikeTimeData;
        private Boolean hasRunned;
        //private KTNeuron[] neuron;
        private INeuron[] neuron;
        private List<ISynapse> synapseList;
        private Int32 nSynapse;
        //private Int32[] indForStimulus;
        private Double stimulusAmp;
        private Double stimulusStdDev;
        private Int32 neuronIndStim;
        private Double poissonRate;
        private ExternalStimulusBase stimulus;
        private Int32 timestep;
        private HomogeneousRand rand;
        private List<Int32> spikeDistribution;
        private List<Int32> timeDistribution;
        private Int32 numOfConnBAGraph;
        private Double rewireProbWSGraph;
        private Boolean hasToAverageInput;
        private Int32 deltaTrainDt;
        private Double dyn_alpha;
        private Double dyn_u;
        private Double dyn_tauJ;
        private Double dyn_dt;
        private KTzV2.Data.NetworkSamplingType samplingType;
        private Double samplingFraction;
        private Int32[] samplingIndices;
        private Int32[] unsampledIndices;
        private Boolean[] countedNeuron;
        private Int32 samplingN;
        private Int32 restIntervals;
        private Int32 timeBin;
        private Int32 NAvalSim;
        private KTzV2.Data.CountVariable countVar;
        private StimulusIndexChoice indChoice;
        public Action<Int32?> RunForAvalanches { get; private set; }
        private Func<Int32> GetNeuronIndexForStim { get; set; }
        private Func<Int32> RunSingleAvalanche { get; set; }
        private Action<Int32> Update { get; set; }
        private Action<Int32> UpdateDuringTransient { get; set; }
        private Action<Double[]> WriteRhoTimeSeriesFile { get; set; }
        public Action<String> WriteDataToFile { get; private set; }
        public Action<String, String> WriteSpikeDistributionToFile { get; private set; }

        private Double[] elAdjMatVal;
        private Func<Int32, Double> GetJValue;
        private Boolean writeRhoTimeSeries; //simulation parameter - writes rho(t) for each J in JRange (only if outAvgMode=OverTime); rho(t) = sum x(t)>0
        private Boolean writeXDataFile;
        private Boolean writeXDataCSVFile;
        private Boolean saveSpikeTimes;
        private KTzV2.Data.Header.OutputFileFormat outputFileFormat;

#if DEBUG
        //private List<Double>[] DEBUG_XData;
        //private List<Int32> DEBUG_TData;
#endif

        public KTzNetworkSimulator(Boolean silent)
        {
            this.rand = new HomogeneousRand();
            Console.WriteLine("Getting parameters...");
            this.getSetOfParameters();
            Console.WriteLine("Setting adjacency matrix...");
            this.createAdjacencyMatrix();

            this.hasRunned = false;

            if (!silent)
            {
                Console.WriteLine("" + "___________");
                Console.WriteLine("");
                Console.WriteLine("" + "    KTz network simulator instanciated");
                Console.WriteLine("");

                Console.WriteLine("" + "-");
            }
        }

        private void createAdjacencyMatrix()
        {
            //Int32 linearSize;
            //if ((networkType == AdjacencyMatrixType.SquareLatticeFreeBC) || (networkType == AdjacencyMatrixType.SquareLatticePeriodicBC))
            //{
            //    linearSize = nNeuronsOnARow;
            //}
            //else
            //{
            //    linearSize = nNeurons;
            //}
            IAdjacencyMatrix<Double> boolMatrixGetter = FGetAdjacencyMatrix.GetMatrixFor(this.networkType, 
                                                                                         this.nNeurons,
                                                                                         new Int32[] { this.Lx, this.Ly, this.Lz },
                                                                                         this.nNeighbours,
                                                                                         this.numOfConnBAGraph,
                                                                                         this.rewireProbWSGraph,
                                                                                         Convert.ToBoolean(KTzHeader.GetPar_Int32(KTzParameters.netDir)), 
                                                                                         KTzHeader.GetPar_String(KTzParameters.netFile)
                                                                                         );
            AMatrix = boolMatrixGetter.BuildAndGetMatrix();
            if (this.networkType == AdjacencyMatrixType.FromFile)
            {
                if (this.AMatrix == null)
                    throw new ArgumentNullException("Could not load and/or create adjacency matrix FromFile");
                else
                {
                    this.nNeurons = this.AMatrix.nRow;
                    this.Lx       = this.nNeurons;
                    this.Ly       = 1;
                    this.Lz       = 1;
                    KTzHeader.SetPar(KTzParameters.Lx, this.Lx);
                    KTzHeader.SetPar(KTzParameters.Ly, this.Ly);
                    KTzHeader.SetPar(KTzParameters.Lz, this.Lz);
                }
            }
            this.GetSampParAndSetSamplingIndices();
        }

        private void getSetOfParameters()
        {
            this.networkType       = (AdjacencyMatrixType)KTzHeader.GetPar_Int32(KTzParameters.netType);
            this.synapseType       = (SynapseType)KTzHeader.GetPar_Int32(KTzParameters.sType);
            this.neuronType        = (NeuronType)KTzHeader.GetPar_Int32(KTzParameters.neuron);
            this.icType            = (InitialConditionType)KTzHeader.GetPar_Int32(KTzParameters.iCond);
            this.countVar          = (KTzV2.Data.CountVariable)KTzHeader.GetPar_Int32(KTzParameters.cVar);
            this.stimulusType      = (StimulusType)KTzHeader.GetPar_Int32(KTzParameters.stimType);
            this.noiseType         = (KTzV2.Synapses.NoiseType)KTzHeader.GetPar_Int32(KTzParameters.noiseType);
            this.couplingParamType = (CouplingParam)KTzHeader.GetPar_Int32(KTzParameters.coupParam);
            this.indChoice         = (StimulusIndexChoice)KTzHeader.GetPar_Int32(KTzParameters.indChoice);
            this.saveSpikeTimes    = KTzHeader.GetPar_Int32(KTzParameters.saveSpikeTimes) == 1;

            /***
             **
             **
             ** Network parameters
             **
             **
             ***/
            this.Lx                = KTzHeader.GetPar_Int32(KTzParameters.Lx);
            this.Ly                = KTzHeader.GetPar_Int32(KTzParameters.Ly);
            this.Lz                = KTzHeader.GetPar_Int32(KTzParameters.Lz);
            this.nNeurons          = this.Lx * this.Ly * this.Lz;//Convert.ToInt32(Math.Pow(Lx, netDim));
            this.netDim            = KTzHeader.GetPar_Int32(KTzParameters.dim);
            this.nNeighbours       = KTzHeader.GetPar_Int32(KTzParameters.nNeigh);
            this.numOfConnBAGraph  = KTzHeader.GetPar_Int32(KTzParameters.nConn);
            this.rewireProbWSGraph = KTzHeader.GetPar_Double(KTzParameters.rewP);

            /***
             **
             **
             ** Simulation parameters
             **
             **
             ***/
            // adjusting RunForAvalanche
            this.ChooseRunForAvalancheMethods();

            // adjusting output writers
            this.ChooseOutputMethods();
            
            // adjusting time step update function
            this.ChooseTimestepUpdateMethod();

            // adjusting constant parameters
            var simTimeScheme   = (SimulationTimeScheme)KTzHeader.GetPar_Int32(KTzParameters.simTimeScheme);
            this.nSteps         = KTzHeader.GetPar_Int32(KTzParameters.nSteps);
            this.nStartStep     = KTzHeader.GetPar_Int32(KTzParameters.nStart);
            this.poissonRate    = KTzHeader.GetPar_Double(KTzParameters.r);
            this.tForStimulus   = KTzHeader.GetPar_Int32(KTzParameters.sStim);
            this.stimulusAmp    = KTzHeader.GetPar_Double(KTzParameters.I);
            this.stimulusStdDev = KTzHeader.GetPar_Double(KTzParameters.IStdDev);
            this.deltaTrainDt   = KTzHeader.GetPar_Int32(KTzParameters.deltaT);
            this.restIntervals  = KTzHeader.GetPar_Int32(KTzParameters.rest);
            this.timeBin        = KTzHeader.GetPar_Int32(KTzParameters.tBin);
            this.NAvalSim       = KTzHeader.GetPar_Int32(KTzParameters.nSim);
            this.neuronIndStim  = KTzHeader.GetPar_Int32(KTzParameters.iStim);

            // adjusting total simulation time if needed
            if ((simTimeScheme == SimulationTimeScheme.ProportionalToPoissonRate) && (this.stimulusType == StimulusType.PoissonProcess))
            {
                var T = this.nStartStep + (int)(10.0D / (double)(this.poissonRate * this.nNeurons));
                if (T > this.nSteps)
                {
                    this.nSteps = T;
                    KTzHeader.SetPar(KTzParameters.nSteps, this.nSteps);
                    Console.WriteLine("WARNING ::: Auto-adjusting total simulation time to nSteps = nStartStep + {0:D}",T);
                }
            }

            /***
             **
             **
             ** Synapse parameters
             **
             **
             ***/
            if (this.couplingParamType == CouplingParam.Homogeneous)
                this.GetJValue = this.GetJValueHomogeneous;
            else
                this.GetJValue = this.GetJValueAdjMatrix;
            this.J           = KTzHeader.GetPar_Double(KTzParameters.J);
            this.tauf        = KTzHeader.GetPar_Double(KTzParameters.tauf);
            this.taug        = KTzHeader.GetPar_Double(KTzParameters.taug);
            this.noiseAmp    = KTzHeader.GetPar_Double(KTzParameters.R);
            this.noiseRatio  = KTzHeader.GetPar_Double(KTzParameters.noiseRatio);
            this.dyn_alpha   = KTzHeader.GetPar_Double(KTzParameters.alpha);
            this.dyn_u       = KTzHeader.GetPar_Double(KTzParameters.u);
            this.dyn_tauJ    = KTzHeader.GetPar_Double(KTzParameters.tauJ);
            this.dyn_dt      = KTzHeader.GetPar_Double(KTzParameters.dt);

            this.hasToAverageInput = KTzHeader.GetPar_Int32(KTzParameters.avgInp) == 1;
            if ((neuronType == NeuronType.KTzLogMF) || (neuronType == NeuronType.KTzMF))
            {
                this.hasToAverageInput = true;
            }

            /***
             **
             **
             ** Neuron parameters
             **
             **
             ***/
            this.K     = KTzHeader.GetPar_Double(KTzParameters.K);
            this.T     = KTzHeader.GetPar_Double(KTzParameters.T);
            this.d     = KTzHeader.GetPar_Double(KTzParameters.d);
            this.l     = KTzHeader.GetPar_Double(KTzParameters.l);
            this.xR    = KTzHeader.GetPar_Double(KTzParameters.xR);
            this.H     = KTzHeader.GetPar_Double(KTzParameters.H);
            this.Q     = KTzHeader.GetPar_Double(KTzParameters.Q);
            this.x0i   = KTzHeader.GetPar_Double(KTzParameters.x0);
            this.y0i   = KTzHeader.GetPar_Double(KTzParameters.y0);
            this.z0i   = KTzHeader.GetPar_Double(KTzParameters.z0);
            this.Theta = KTzHeader.GetPar_Double(KTzParameters.Theta);


            this.doWriteDif = false;
            if (KTzHeader.GetPar_Int32(KTzParameters.wDif) == 1)
                Console.WriteLine("WARNING: parameter wDif has no effect anymore");
        }

        private void ChooseRunForAvalancheMethods()
        {
            var dynType            = (KTzV2.DynamicsSimType)KTzHeader.GetPar_Int32(KTzParameters.dynType);
            var simType            = (KTzV2.SimulationType)KTzHeader.GetPar_Int32(KTzParameters.simType);

            // if we want to run the dynamics (instead of bifurcation)
            // and the stimulus IS NOT DeltaWhenInactive
            // we force the simulation to be "ContinuousTime"
            if ((this.stimulusType != StimulusType.DeltaWhenInactive) && (simType == KTzV2.SimulationType.Dynamics))
            {
                dynType = KTzV2.DynamicsSimType.ContinuousTime;
                KTzHeader.SetPar(KTzParameters.dynType, (Int32)KTzV2.DynamicsSimType.ContinuousTime);
                Console.WriteLine("WARNING: Forcing dynType == ContinuousTime because stimType != DeltaWhenInactive");
            }

            if (dynType == KTzV2.DynamicsSimType.NetworkReset)
            {
                this.stimulusType = StimulusType.DeltaWhenInactive;
                KTzHeader.SetPar(KTzParameters.stimType, (Int32)StimulusType.DeltaWhenInactive);
                Console.WriteLine("WARNING: Forcing stimType == DeltaWhenInactive because dynType == NetworkReset");
            }

            if (dynType == KTzV2.DynamicsSimType.ContinuousTime)
            {
                this.RunForAvalanches = RunForSpikeDistribution;
            }
            else if (dynType == KTzV2.DynamicsSimType.NetworkReset)
            {
                this.RunForAvalanches = RunForSpikeDistReset;
                this.stimulusType     = StimulusType.Delta;
                KTzHeader.SetPar(KTzParameters.stimType, (int)this.stimulusType);
            }

            // adjusting RunSingleAvalanche
            if (this.countVar == CountVariable.NumberOfNeurons)
                this.RunSingleAvalanche = this.RunSingleAvalancheForNeurons;
            else if (this.countVar == CountVariable.NumberOfSpikes)
                this.RunSingleAvalanche = this.RunSingleAvalancheForSpike;
        }

        private void ChooseOutputMethods()
        {
            this.outputFileFormat   = (KTzV2.Data.Header.OutputFileFormat)KTzHeader.GetPar_Int32(KTzParameters.oFileFormat);
            if (this.outputFileFormat == OutputFileFormat.txt)
            {
                this.WriteRhoTimeSeriesFile       = this.WriteRhoTimeSeriesFileTxt;
                this.WriteDataToFile              = this.WriteDataToFileTxt;
                this.WriteSpikeDistributionToFile = this.WriteSpikeDistributionToFileTxt;
            }
            else
            {
                this.WriteRhoTimeSeriesFile       = this.WriteRhoTimeSeriesFileMat;
                this.WriteDataToFile              = this.WriteDataToFileMat;
                this.WriteSpikeDistributionToFile = this.WriteSpikeDistributionToFileMat;
            }
            this.writeRhoTimeSeries = KTzHeader.GetPar_Int32(KTzParameters.writeRhoTS) == 1;
            this.writeXDataFile     = KTzHeader.GetPar_Int32(KTzParameters.wData) == 1;
            this.writeXDataCSVFile  = KTzHeader.GetPar_Int32(KTzParameters.wCSV) == 1;
        }

        private void ChooseTimestepUpdateMethod()
        {
            if (this.saveSpikeTimes)
            {
                if (this.indChoice == StimulusIndexChoice.MultipleRandom)
                {
                    this.Update                = this.UpdateForMultipleStimSaveSpike;
                    this.UpdateDuringTransient = this.UpdateForMultipleStim;
                }
                else
                {
                    this.Update                = this.UpdateSingleStimSaveSpike;
                    this.UpdateDuringTransient = this.UpdateSingleStim;
                }
            }
            else
            {
                if (this.indChoice == StimulusIndexChoice.MultipleRandom)
                {
                    this.Update                = this.UpdateForMultipleStim;
                    this.UpdateDuringTransient = this.UpdateForMultipleStim;
                }
                else
                {
                    this.Update                = this.UpdateSingleStim;
                    this.UpdateDuringTransient = this.UpdateSingleStim;
                }
            }
        }

        private Int32 GetNeuronIndForStimAll()
        {
            return -1;
        }

        private Int32 GetNeuronIndForStimRandom()
        {
            return (Int32)(this.rand.GetRandomLT1() * this.nNeurons);
        }

        private Int32 GetNeuronIndForStimSimple()
        {
            return this.neuronIndStim;
        }

        private void GetSampParAndSetSamplingIndices()
        {
            this.samplingFraction = KTzHeader.GetPar_Double(KTzParameters.sampFrac);
            this.samplingType     = (KTzV2.Data.NetworkSamplingType)KTzHeader.GetPar_Int32(KTzParameters.samp);
            this.samplingN        = (Int32)(this.samplingFraction * this.nNeurons);
            if (this.samplingN == this.nNeurons)
                this.samplingType = NetworkSamplingType.Full;

            Int32 i;
            if (this.samplingType == KTzV2.Data.NetworkSamplingType.Full)
            {
                this.unsampledIndices = new Int32[] {};
                this.samplingN        = this.nNeurons;
                this.samplingIndices  = new Int32[this.samplingN];
                i = 0;
                while (i < this.samplingN)
                {
                    this.samplingIndices[i] = i;
                    i++;
                }
            }
            else if (this.samplingType == KTzV2.Data.NetworkSamplingType.Partial)
            {
                this.samplingN       = (Int32)(this.samplingFraction * this.nNeurons);
                this.samplingIndices = new Int32[this.samplingN];
                List<Int32> selected = new List<Int32>();
                Int32 n;
                i = 0;
                while (i < this.samplingN)
                {
                    do
                    {
                        n = (Int32)(this.rand.GetRandomLT1() * this.nNeurons); // choose a random index different of the already taken
                    }
                    while (selected.Contains(n));
                    this.samplingIndices[i] = n;
                    selected.Add(n);
                    i++;
                }
                selected = null;
                this.unsampledIndices = new Int32[this.nNeurons - this.samplingN];
                var diff = Enumerable.Range(0, this.nNeurons).Except(this.samplingIndices);
                int k = 0;
                foreach (int d in diff)
                {
                    this.unsampledIndices[k] = d;
                    k++;
                }
            }
            else
            {
                throw new ArgumentOutOfRangeException("Unrecognized sampling type! " + this.samplingType.ToString());
            }
        }

        private void prepareNeurons()
        {
            Int32 i, j, k;
            //Console.WriteLine(lb + " * Generating initial conditions...");
            this.x0p = new Double[nNeurons];
            this.y0p = new Double[nNeurons];
            this.z0p = new Double[nNeurons];
            this.icGetter = new GetICInstance();
            i = 0;
            try
            {
                this.icGetter = FGetInitCond.use(this.icType, 1, nNeurons);
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
                Console.WriteLine("Using type 3 as default initial condition!");
                this.icType = (InitialConditionType)3;
                this.icGetter = FGetInitCond.use(this.icType, 1, nNeurons);
            }
            while (i < this.nNeurons)
            {
                if (this.icType == InitialConditionType.FromXMLFile)
                {
                    this.x0p[i] = icGetter.getX(0, i);
                    this.y0p[i] = icGetter.getY(0, i);
                    this.z0p[i] = icGetter.getZ(0, i);
                }
                else
                {
                    this.x0p[i] = icGetter.getX(x0i);
                    this.y0p[i] = icGetter.getY(y0i);
                    this.z0p[i] = icGetter.getZ(z0i);
                }
                i++;
            }

            this.countedNeuron = new Boolean[this.nNeurons];
            this.neuron = new INeuron[nNeurons]; // array of neurons
            i = 0;
            while (i < this.nNeurons)
            {
                //neuron[i] = new KTzNeuron(i, K, T, d, l, xR, x0[i], y0[i], z0[i]);
                // I have to get the parameters here again because these parameters may have quenched disorder
                this.K         = KTzHeader.GetPar_Double(KTzParameters.K);
                this.T         = KTzHeader.GetPar_Double(KTzParameters.T);
                this.d         = KTzHeader.GetPar_Double(KTzParameters.d);
                this.l         = KTzHeader.GetPar_Double(KTzParameters.l);
                this.xR        = KTzHeader.GetPar_Double(KTzParameters.xR);
                this.H         = KTzHeader.GetPar_Double(KTzParameters.H);
                this.Q         = KTzHeader.GetPar_Double(KTzParameters.Q);
                this.Theta     = KTzHeader.GetPar_Double(KTzParameters.Theta);
                this.neuron[i] = NeuronFactory.GetKTNeuron(this.neuronType,
                                                           new NeuronParam(i, this.K, this.T, this.d, this.l, this.xR, this.H, this.x0p[i], this.y0p[i], this.z0p[i], this.Theta, Q: this.Q)
                                                           );
                this.countedNeuron[i] = false;
                i++;
            }

            this.synapseList = new List<ISynapse>();
            Int32[] preNeuronsInd;

            // adding the connections to the neurons
            UInt64 seed = 0;
            i = 0;
            while (i < this.nNeurons)
            {
                // getting the indeces of the pre-synaptic neurons of the i-th neuron
                preNeuronsInd = this.AMatrix.GetNonEmptyRowsInd(i, out this.elAdjMatVal);
                j = 0;
                k = preNeuronsInd.Length; // amount of neurons connected to this one
                while (j < k)
                {
                    // I have to retrieve these parameters because they may be subject to quenched disorder
                    this.J          = this.GetJValue(j);
                    this.tauf       = KTzHeader.GetPar_Double(KTzParameters.tauf);
                    this.taug       = KTzHeader.GetPar_Double(KTzParameters.taug);
                    this.noiseAmp   = KTzHeader.GetPar_Double(KTzParameters.R);
                    this.noiseRatio = KTzHeader.GetPar_Double(KTzParameters.noiseRatio);
                    this.dyn_alpha  = KTzHeader.GetPar_Double(KTzParameters.alpha);
                    this.dyn_u      = KTzHeader.GetPar_Double(KTzParameters.u);
                    this.dyn_tauJ   = KTzHeader.GetPar_Double(KTzParameters.tauJ);
                    this.synapseList.Add(
                            SynapseFactory.GetSynapse(
                                this.synapseType, new SynapseParam(this.neuron[preNeuronsInd[j]], this.neuron[i],
                                                                   this.tauf, this.taug, this.J, this.noiseAmp, this.noiseRatio,
                                                                   seed + (UInt64)DateTime.Now.Ticks,
                                                                   this.noiseType,
                                                                   (this.dyn_alpha / this.dyn_u), this.dyn_tauJ, this.dyn_alpha, this.dyn_u, this.dyn_dt)
                            )
                        );
                    this.neuron[i].AddInput(this.synapseList.Last());
                    seed++;
                    j++;
                }
                i++;
            }
            this.nSynapse = this.synapseList.Count;
        }

        private Double GetJValueHomogeneous(Int32 i)
        {
            return KTzHeader.GetPar_Double(KTzParameters.J);
        }

        private Double GetJValueAdjMatrix(Int32 i)
        {
            return this.elAdjMatVal[i];
        }

        private void prepareStimulus()
        {
            // neurons are already prepared at this point
            if (this.neuronIndStim == -1)
            {
                if (this.indChoice == StimulusIndexChoice.SquareCenter)
                {
                    this.neuronIndStim = this.GetSquareCenterIndex();
                    this.GetNeuronIndexForStim = GetNeuronIndForStimSimple;
                }
                else if (this.indChoice == StimulusIndexChoice.Fixed)
                {
                    this.neuronIndStim = this.GetSquareCenterIndex();
                    this.GetNeuronIndexForStim = GetNeuronIndForStimSimple;
                }
                else if (this.indChoice == StimulusIndexChoice.MostConnected)
                {
                    this.neuronIndStim = this.GetMostConnectedNeuronIndex();
                    this.GetNeuronIndexForStim = GetNeuronIndForStimSimple;
                }
                else if (this.indChoice == StimulusIndexChoice.Random)
                {
                    this.GetNeuronIndexForStim = GetNeuronIndForStimRandom;
                }
                else if (this.indChoice == StimulusIndexChoice.MultipleRandom)
                {
                    this.GetNeuronIndexForStim = GetNeuronIndForStimAll;
                }
                else
                {
                    throw new ArgumentOutOfRangeException("Unknown StimulusIndexChoice");
                }
            }
            else
            {
                KTzHeader.SetPar(KTzParameters.indChoice, (Int32)StimulusIndexChoice.Fixed);
                this.GetNeuronIndexForStim = GetNeuronIndForStimSimple;
            }

            // each time we simulate, we create a new stimulus vector
            this.stimulus = FExternalStimulus.Get(this.stimulusType, this.tForStimulus, this.stimulusAmp, this.deltaTrainDt, this.nSteps / this.deltaTrainDt, 1.0 / this.poissonRate, this.stimulusStdDev);
        }

        private Int32 GetSquareCenterIndex()
        {
            return (this.networkType == AdjacencyMatrixType.SquareLatticeFreeBC || this.networkType == AdjacencyMatrixType.SquareLatticePeriodicBC) ? ((Int32)Math.Floor(this.Lx / 2.0D) * this.Ly + (Int32)Math.Floor(this.Lx / 2.0D)) : (this.nNeurons / 2);
        }

        private Int32 GetMostConnectedNeuronIndex()
        {
            return this.neuron.Aggregate(this.neuron[0], (curMax, next) => next.Input.Count > curMax.Input.Count ? next : curMax).Index;
        }

        public void SetParam(KTzParameters param, Double value)
        {
            if (!KTzHeader.HasParam(param))
            {
                throw new Exception("Specified parameter does not exist!");
            }
            KTzHeader.SetPar(param, value);
            hasRunned = false;
            getSetOfParameters();
            if (param == KTzParameters.rewP)
            {
                Console.WriteLine("Resetting adjacency matrix...");
                createAdjacencyMatrix();
            }
        }

        public void SetParam(KTzParameters param, Int32 value)
        {
            if (!KTzHeader.HasParam(param))
            {
                throw new Exception("Specified parameter does not exist!");
            }
            KTzHeader.SetPar(param, value);
            hasRunned = false;
            getSetOfParameters();
            if ((param == KTzParameters.Lx) ||
                (param == KTzParameters.Ly) ||
                (param == KTzParameters.Lz) ||
                (param == KTzParameters.netType) ||
                (param == KTzParameters.nNeigh) ||
                (param == KTzParameters.nConn) ||
                (param == KTzParameters.dim))
            {
                Console.WriteLine("Resetting adjacency matrix...");
                createAdjacencyMatrix();
            }
        }

        public void ClearSimulation()
        {
            this.xData = null;
            //this.tData = null;
            this.hasRunned = false;

            if (icType != InitialConditionType.FromLastSimulation)
            {
                this.ResetNeuronsAndSynapses();
            }
        }

        private void ResetNeuronsAndSynapses()
        {
            this.timestep = -1;
            int i = 0;
            while (i < this.nNeurons)
            {
                this.neuron[i].ResetIC(new Double[3] { this.x0p[i], this.y0p[i], this.z0p[i] });
                i++;
            }
            this.synapseList.ForEach(ss => ss.ResetSignal());
        }

        private void SaveSpikeTime(int neuronInd, int time)
        {
            this.spikeNeuronData.Add(neuronInd);
            this.spikeTimeData.Add(time);
        }

        private void UpdateSingleStimSaveSpike(Int32 i = 0)
        {
            // evaluating the output synapses of each neuron
            Int32 k = 0, j;
            while (k < this.nSynapse)
            {
                synapseList[k].Evolve();
                k++;
            }
            // evaluating the membrane potential of each neuron
            k = 0;
            while (k < this.samplingIndices.Length)
            {
                j = this.samplingIndices[k];
                if (j == i)
                    neuron[j].Evolve();
                else
                    neuron[j].Evolve(stimulus.GetSignal(this.timestep));
                if (neuron[j].SpikeDetector())
                    this.SaveSpikeTime(j, this.timestep);
                k++;
            }
            k = 0;
            while (k < this.unsampledIndices.Length)
            {
                j = this.unsampledIndices[k];
                if (j == i)
                    neuron[j].Evolve();
                else
                    neuron[j].Evolve(stimulus.GetSignal(this.timestep));
                k++;
            }
        }

        /// <summary>
        /// updates this network 1 timestep and the given neuron to stimulate
        /// if synapses have no dynamics at all (Pulse coupling or GapJunction)
        /// </summary>
        /// <param name="i">index of the neuron to stimulate</param>
        private void UpdateSingleStim(Int32 i = 0)
        {
            // evaluating the output synapses of each neuron
            Int32 j = 0;
            while (j < this.nSynapse)
            {
                synapseList[j].Evolve();
                j++;
            }
            // evaluating the membrane potential of each neuron
            j = 0;
            while (j < i)
            {
                neuron[j].Evolve();
                j++;
            }
            neuron[j].Evolve(stimulus.GetSignal(this.timestep));
            j++;
            while (j < this.nNeurons)
            {
                neuron[j].Evolve();
                j++;
            }
        }

        private void UpdateForMultipleStimSaveSpike(Int32 i = 0)
        {
            // evaluating the output synapses of each neuron
            Int32 k = 0, j;
            while (k < this.nSynapse)
            {
                synapseList[k].Evolve();
                k++;
            }
            // evaluating the membrane potential of each neuron
            k = 0;
            while (k < this.samplingIndices.Length)
            {
                j = this.samplingIndices[k];
                neuron[j].Evolve(stimulus.GetSignal(this.timestep));
                if (neuron[j].SpikeDetector())
                    this.SaveSpikeTime(j, this.timestep);
                k++;
            }
            k = 0;
            while (k < this.unsampledIndices.Length)
            {
                neuron[this.unsampledIndices[k]].Evolve(stimulus.GetSignal(this.timestep));
                k++;
            }
        }

        private void UpdateForMultipleStim(Int32 i = 0)
        {
            // evaluating the output synapses of each neuron
            Int32 j = 0;
            while (j < this.nSynapse)
            {
                synapseList[j].Evolve();
                j++;
            }
            // evaluating the membrane potential of each neuron
            j = 0;
            while (j < this.nNeurons)
            {
                neuron[j].Evolve(stimulus.GetSignal(this.timestep));
                j++;
            }
        }

        /// <summary>
        /// prepares this simulation... prepares the data vars, the neurons, the connections and the stimulus
        /// </summary>
        private void prepareToRunWithData(bool recordXData = true)
        {
            this.prepareToRunWithData(this.nSteps - this.nStartStep + 1, recordXData);
        }

        /// <summary>
        /// prepares the simulation to store an defined amount of data
        /// </summary>
        /// <param name="amountOfData">the amount of data to be stored (amount of timesteps that simulation will run)</param>
        private void prepareToRunWithData(Int32 amountOfData, bool recordXData = true)
        {
            this.prepareToRun();
            if (recordXData)
            {
                this.xData = new Double[this.nNeurons][];
                //this.tData = new Int32[amountOfData];
                Int32 i = 0;
                while (i < this.nNeurons)
                {
                    this.xData[i] = new Double[amountOfData];
                    i++;
                }
            }

            this.spikeNeuronData   = new List<int>();
            this.spikeTimeData     = new List<int>();
            this.spikeDistribution = new List<int>();
            this.timeDistribution  = new List<int>();
        }

        private void prepareToRun()
        {
            prepareNeurons();

            prepareStimulus();
        }

        /// <summary>
        /// runs the simulation for a random neuron index for stimulus
        /// </summary>
        public void Run(Int32? indToStimulate = null, bool recordXData = true)
        {
            if (hasRunned)
            {
                throw new Exception("Simulation already runned!");
            }

            this.prepareToRunWithData(recordXData);

            if (indToStimulate.HasValue)
            {
                KTzHeader.SetPar(KTzParameters.indChoice, (Int32)StimulusIndexChoice.Fixed);
                this.neuronIndStim = indToStimulate.Value;
            }

            this.timestep = -1;
            while (this.timestep < this.nStartStep)
            {
                this.timestep++;
                this.UpdateDuringTransient(GetNeuronIndexForStim());
            }
            if (recordXData)
            {
                if (this.stimulusType == StimulusType.DeltaWhenInactive)
                    this.RunDeltaInactive(this.nSteps - this.nStartStep + 1);
                else
                    this.RunAndRecordTimeLength(this.nSteps - this.nStartStep + 1);
            }
            else
                this.RunTimeLength(this.nSteps - this.nStartStep + 1);

            hasRunned = true;
        }

        private void RunDeltaInactive(Int32 timeLength)
        {
            Double[][] xDataTemp = new Double[this.nNeurons][];
            this.prepareToRunWithData(this.timeBin);
            Int32 t = 0, countTemp = 0;
            Int32 tMax = this.timeBin * (timeLength / this.timeBin); // total time
            for (int i = 0; i < this.nNeurons; i++)
            {
                xDataTemp[i] = new Double[tMax];
            }

            this.timestep = -1;
            while (t < tMax)
            {
                // running the simulation on this bin
                this.RunAndRecordTimeLength(this.timeBin, ref xDataTemp);

                countTemp = this.CountNonConsecutiveDataPeaks(0, this.timeBin); // counts spikes only in the sampled data

                // checks the full data for spikes - needed for deciding whether to stimulate the network
                if ((countTemp == 0) && this.IsThereNoActivityInUnsampledData()) // if there's no activity in complete sampled data
                {
                    if (this.stimulus.timeForStimulus < t)
                        this.stimulus.SetNextStimulusTime(1+this.timestep + this.restIntervals * this.timeBin); // set the next stimulus timestep
                }

                // setting the start of the next bin
                t += this.timeBin;
            }
            this.xData = xDataTemp;
        }

        /// <summary>
        /// runs this simulation in order to get the spike count and returns it, storing it
        /// </summary>
        /// <param name="tWindow">the time bin length to count spikes</param>
        public ThermoStatistics RunForSpikeCountTemporalDynamics(Int32 indexForStimulus = -1, Int32 tWindow = 1000)
        {
            Int32? indSt;
            if (indexForStimulus > -1)
                indSt = indexForStimulus;
            else
                indSt = null;

            //
            //
            // simulation has to run many times
            //
            //
            if (this.saveSpikeTimes)
                this.Run(indSt, recordXData: false);
            else
                this.Run(indSt);

            Double[] rho_t;
            Double rho, rho2, rho4, rhoMax;
            Int32 t_total;
            if (this.saveSpikeTimes)
                this.GenerateRhoTimeSeriesFromSpikeData(out rho_t, out rho, out rho2, out rho4, out rhoMax, out t_total);
            else
                this.GenerateRhoTimeSeriesFromXData(out rho_t, out rho, out rho2, out rho4, out rhoMax, out t_total);


            if (this.writeRhoTimeSeries)
                this.WriteRhoTimeSeriesFile(rho_t);

            if (this.writeXDataFile)
                this.WriteDataToFile("");

            //
            //
            // mean has to be taken over many realizations
            //
            //
            return new ThermoStatistics(rho, rho2, rho4, rho_t, rhoMax, t_total);
        }

        private void GenerateRhoTimeSeriesFromSpikeData(out Double[] rho_t, out Double rho, out Double rho2, out Double rho4, out Double rhoMax, out Int32 t_total)
        {
            Double n_sample = (Double)this.samplingIndices.Length;
            t_total         = -1;
            Double rho_tmp  = 0.0D;
            rho             = 0.0D;
            rho2            = 0.0D;
            rho4            = 0.0D;
            rhoMax          = 0.0D;
            var rho_t_list  = new List<Double>();
            var n           = 1; // number of neurons that spiked at t0
            var t0          = this.spikeTimeData[0]; // instant of a spike
            var t0_prev     = -1;
            for (var k = 1; k < this.spikeTimeData.Count; k++)
            {
                // if there is some spike at t0
                if (this.spikeTimeData[k] == t0)
                {
                    n++; // then another neuron spiked at this time
                }
                else // otherwise, we restart counting the spike for a new t0
                {
                    // calculating the fraction of spikes in time t0
                    rho_tmp = (Double)n / n_sample;

                    // resetting spike time
                    t0_prev = t0;
                    n       = 1; // there is already one spike at this new t0
                    t0      = this.spikeTimeData[k]; // this is the new spike instant

                    // recording the fraction of spikes in the sampled data
                    rho_t_list.Add(rho_tmp);
                    // filling with zeros if needed
                    if ((t0 - t0_prev) > 1) // if the next spike did not happen immediately
                        rho_t_list.AddRange(Enumerable.Range(t0_prev,t0-1-t0_prev).Select(k=>0.0D)); // we add a bunch of zeros since the network remained silent during (t0-1-t0_prev) time steps
                    
                    // calculating means
                    // the zeros added to the list don't change the sums below
                    rho  += rho_tmp;
                    rho2 += rho_tmp * rho_tmp;
                    rho4 += rho_tmp * rho_tmp * rho_tmp * rho_tmp;
                }
            }
            rho_t   = rho_t_list.ToArray();
            t_total = rho_t.Length;
            rho    /= (Double)t_total;
            rho2   /= (Double)t_total;
            rho4   /= (Double)t_total;
        }

        private void GenerateRhoTimeSeriesFromXData(out Double[] rho_t, out Double rho, out Double rho2, out Double rho4, out Double rhoMax, out Int32 tt)
        {
            //Func<Double, Double> linearTransf;
            //if (this.neuronType.ToString().Contains("KT"))
            //    linearTransf = x => (x + 1.0D) / 2.0D;
            //else
            //    linearTransf = x => x;

            Int32 t = 1, i, k;
            Int32 n_sample = this.samplingIndices.Length;
            tt             = this.xData[0].Length;
            rho            = 0.0D;
            rho2           = 0.0D;
            rho4           = 0.0D;
            rhoMax         = 0.0D;
            rho_t          = new Double[tt];
            while (t < tt)
            {
                rho_t[t] = 0.0D;
                i = 0;
                while (i < n_sample) //this.nNeurons)
                {
                    k = this.samplingIndices[i];
                    // this was already commented //rho_t[t] += linearTransf(this.xData[i][t]); // transforming xData from [-1;1] to [0;1] if needed
                    //rho_t[t] += (this.xData[i][t] > 0.0D ? 1.0D : 0.0D); // summing up only the spikes
                    rho_t[t] += this.neuron[k].SpikeDetector(this.xData[k][t],this.xData[k][t-1]) ? 1.0D : 0.0D; // summing up only the spikes
                    i++;
                }
                rho_t[t] /= n_sample;
                rho      += rho_t[t];
                rho2     += rho_t[t] * rho_t[t];
                rho4     += rho_t[t] * rho_t[t] * rho_t[t] * rho_t[t];
                if (rho_t[t] > rhoMax)
                    rhoMax = rho_t[t];
                t++;
            }
            rho  /= (Double)tt;
            rho2 /= (Double)tt;
            rho4 /= (Double)tt;
        }

        private void WriteRhoTimeSeriesFileMat(Double[] s)
        {
            String fileName = KTzHeader.GetPar_String(KTzParameters.oFile);
            fileName = KTzHeader.CheckAndGetFileName(fileName + (fileName.Length > 0 ? "_" : "") + "rho_timeseries.mat");
            var matDataBuilder = new MatFileHandler.DataBuilder();
            Dictionary<KTzParamGroups, MatFileHandler.IVariable> simParams = KTzHeader.GetAllParamPairsAsMatlabStruct(matDataBuilder);
            using (var fileStream = new System.IO.FileStream(fileName, System.IO.FileMode.Create))
            {
                Console.WriteLine("" + " * Writing to file: {0}", fileName);
                var writer = new MatFileHandler.MatFileWriter(fileStream);
                writer.Write(
                    matDataBuilder.NewFile(simParams.Values.Concat(new[]{
                         matDataBuilder.NewVariable("rho", matDataBuilder.NewArray(s, s.Length, 1)),
                         matDataBuilder.NewVariable("time", matDataBuilder.NewArray(Enumerable.Range(0, s.Length).ToArray(), s.Length, 1)),
                         matDataBuilder.NewVariable("file_info", matDataBuilder.NewCharArray("rho: spike count time series"))
                    }))
                );
            }
        }/**/

        private void WriteRhoTimeSeriesFileTxt(Double[] s)
        {
            // adjusting the output file header
            String fHeader = "# Density of active sites (rho) versus time" + Environment.NewLine;
            fHeader += "#-" + Environment.NewLine;
            fHeader += this.GetOutputFileHeader();
            fHeader += "#-" + Environment.NewLine;
            String[] oFileColNames = new String[] { "time", "rho" };

            String oFileName = KTzHeader.GetPar_String(KTzParameters.oFile);
            oFileName += (oFileName.Length > 0 ? "_" : "") + "rho_timeseries.dat";
            Console.WriteLine("" + " * Writing to file: {0}", oFileName);
            System.IO.StreamWriter sw = KTzHeader.CreateOutTxtFile(oFileName, fHeader, oFileColNames);
            int i = 0;
            while (i < s.Length)
            {
                sw.WriteLine("{0}\t{1:0.00000000e+000}", i, s[i]);
                i++;
            }
            KTzHeader.CloseOutTxtFile(sw);
        }

        /// <summary>
        /// runs this simulation NAvalSim times in order to get the spike count and returns it, storing it
        /// </summary>
        /// <param name="tWindow">the time bin length to count spikes</param>
        public ThermoStatistics RunNTimesForSpikeCount(Int32 indexForStimulus = -1, Int32 tWindow = 1000)
        {
            if (tWindow <= 0)
            {
                throw new ArgumentException("invalid binLength");
            }

            // preparing to run
            this.prepareToRunWithData(tWindow);

            if (indexForStimulus != -1)
                this.neuronIndStim = indexForStimulus;

            // t = beginning of time bin, tNext = end of time bin, count = current count of spikes in consecutive time bins, countTemp = count of spikes in the current time bin
            Double c = 0.0D, c2 = 0.0D, c4 = 0.0D, cMax = 0.0D;
            Double[] cPerSim = new Double[this.NAvalSim];
            Int32 tMax = tWindow * (this.nSteps / tWindow); // total time
            Int32 t;//, tNext;//, countTemp = 0;
            Int32 l = 0;

            while (l < this.NAvalSim)
            {
                cPerSim[l] = 0.0D;

                this.ResetNeuronsAndSynapses();
                t = 0;
                while (t < tMax)
                {
                    this.RunAndRecordTimeLength(tWindow);

                    // counting the spikes for this bin
                    cPerSim[l] += (Double)this.CountNonConsecutiveDataPeaks(0, tWindow);

                    // setting the start of the next bin
                    t = t + tWindow;
                }

                if (cPerSim[l] > cMax)
                    cMax = cPerSim[l];
                c += cPerSim[l];
                c2 += cPerSim[l] * cPerSim[l];
                c4 += cPerSim[l] * cPerSim[l] * cPerSim[l] * cPerSim[l];

                l++;
            }
            return new ThermoStatistics(c / (Double)this.NAvalSim, c2 / (Double)this.NAvalSim, c4 / (Double)this.NAvalSim, cPerSim, cMax, this.NAvalSim);
        }

        /// <summary>
        /// runs this simulation in order to get the spike count and returns it, storing it
        /// </summary>
        /// <param name="tWindow">the time bin length to count spikes</param>
        /*public Int32 RunForNeuronCount(Int32 indexForStimulus)
        {
            if (indexForStimulus != -1)
            {
                this.Run(indexForStimulus);
            }
            else
            {
                this.Run();
            }
            return this.CountSpikingNeurons();
        }*/

        public ThermoStatistics RunNAvalanches()
        {
            if (this.countVar == KTzV2.Data.CountVariable.NumberOfSpikes)
            {
                this.prepareToRunWithData(this.timeBin);
            }
            else if (this.countVar == KTzV2.Data.CountVariable.NumberOfNeurons)
            {
                this.nStartStep = 0;
                this.prepareToRunWithData(this.nSteps);
            }
            Double m = 0, mSD = 0, u = 0, mMax = 0;
            Double[] countPerSim = new Double[this.NAvalSim];
            Int32 l = 0; // performing nSim simulations in order to make get a mean value for spikeCountTot and neuronCount
            while (l < this.NAvalSim)
            {
                // simulating
                this.ResetNeuronsAndSynapses();
                countPerSim[l] = this.RunSingleAvalanche();
                if (countPerSim[l] > mMax)
                    mMax = countPerSim[l];
                m += countPerSim[l];
                mSD += countPerSim[l] * countPerSim[l];
                u += countPerSim[l] * countPerSim[l] * countPerSim[l] * countPerSim[l];

                l++;
            }
            // I can't parallelize this loop because each runs updates the same set of neurons... I would need to create a new set of neurons and synapses for each step of the loop
            //Parallel.For(0, this.NAvalSim, l =>
            //    {
            //        this.ResetNeuronsAndSynapses();
            //        countPerSim[l] = this.RunSingleAvalanche();
            //        if (countPerSim[l] > mMax)
            //            mMax = countPerSim[l];
            //        m += countPerSim[l];
            //        mSD += countPerSim[l] * countPerSim[l];
            //        u += countPerSim[l] * countPerSim[l] * countPerSim[l] * countPerSim[l];
            //    });
            return new ThermoStatistics(m / (Double)this.NAvalSim, mSD / (Double)this.NAvalSim, u / (Double)this.NAvalSim, countPerSim, mMax, this.nNeurons);
        }

        //private Int32 RunSingleAvalancheForNeurons()
        //{
        //    this.RunAndRecordTimeLength(this.nSteps);
        //    return this.CountSpikingNeurons();
        //}
        private Int32 RunSingleAvalancheForNeurons()
        {
            Boolean noActivityInFullBinData;
            Int32 t = 0;
            //Int32 avSize = 0;
            //avDuration = 0;
            Int32 tMax = this.timeBin * (this.nSteps / this.timeBin); // total time
            this.countedNeuron = new Boolean[this.nNeurons];

            this.timestep = -1;
            while (t < tMax)
            {
                // running the simulation on this bin
                this.RunAndRecordTimeLength(this.timeBin);

                this.ChangeNeuronCountedState();

                noActivityInFullBinData = (this.CountNonConsecutiveDataPeaks(0, this.timeBin) == 0) && this.IsThereNoActivityInUnsampledData(); // checks the full data for spikes - needed for deciding whether to stimulate the network
                if (noActivityInFullBinData) // if there's no activity in complete sampled data
                {
                    break;
                }

                // setting the start of the next bin
                t += this.timeBin;
            }
            return this.countedNeuron.Where((e) => e).Count();
        }

        private Int32 RunSingleAvalancheForSpike()
        {
            Boolean noActivityInFullBinData;
            Int32 t = 0, countTemp = 0;
            Int32 avSize = 0;
            //avDuration = 0;
            Int32 tMax = this.timeBin * (this.nSteps / this.timeBin); // total time

            this.timestep = -1;
            while (t < tMax)
            {
                // running the simulation on this bin
                this.RunAndRecordTimeLength(this.timeBin);

                countTemp = this.CountNonConsecutiveDataPeaks(0, this.timeBin); // counts spikes only in the sampled data
                avSize += countTemp;

                noActivityInFullBinData = (countTemp == 0) && this.IsThereNoActivityInUnsampledData(); // checks the full data for spikes - needed for deciding whether to stimulate the network
                if (noActivityInFullBinData) // if there's no activity in complete sampled data
                {
                    break;
                }

                // setting the start of the next bin
                t += this.timeBin;
            }
            return avSize;
        }

        private void RunAndRecordNumOfSpikesInTimeLength(Int32 timeLength)
        {
            Int32 tt = 1, i = 0;
            // saving last state in previous window so we don't miss any spike in between time windows
            while (i < this.nNeurons)
            {
                this.xData[i][0] = neuron[i].x;
                i++;
            }
            while (tt < timeLength)
            {
                this.timestep++;
                this.Update(GetNeuronIndexForStim());

                i = 0;
                while (i < this.nNeurons)
                {
                    this.xData[i][tt] = neuron[i].x;
                    i++;
                }
                tt++;
            }
        }

        private void RunAndRecordTimeLength(Int32 timeLength)
        {
            Int32 tt = 1, i = 0;
            // saving last state in previous window so we don't miss any spike in between time windows
            while (i < this.nNeurons)
            {
                this.xData[i][0] = neuron[i].x;
                i++;
            }
            while (tt < timeLength)
            {
                this.timestep++;
                this.Update(GetNeuronIndexForStim());

                i = 0;
                while (i < this.nNeurons)
                {
                    this.xData[i][tt] = neuron[i].x;
                    i++;
                }
                tt++;
            }
        }

        private void RunTimeLength(Int32 timeLength)
        {
            int t = -1;
            while (t < timeLength)
            {
                t++;
                this.timestep++;
                this.Update(GetNeuronIndexForStim());
            }
        }

        private void RunAndRecordTimeLength(Int32 timeLength, ref Double[][] cData)
        {
            Int32 tt = 1, i = 0;
            // saving last state in previous window so we don't miss any spike in between time bins
            if (this.timestep >= 0)
                this.timestep--; // if the simulation already ran for at least 1 time bin, then we retrocede 1 timestep to record again the last state of the previous time window
            while (i < this.nNeurons)
            {
                this.timestep++;
                this.xData[i][0] = neuron[i].x;
                cData[i][this.timestep] = neuron[i].x;
                i++;
            }
            while (tt < timeLength)
            {
                this.timestep++;
                this.Update(GetNeuronIndexForStim());

                i = 0;
                while (i < this.nNeurons)
                {
                    this.xData[i][tt] = neuron[i].x;
                    cData[i][this.timestep] = neuron[i].x;
                    i++;
                }
                tt++;
            }
        }

        private void RunForSpikeDistReset(Int32? tBinInput = 14)
        {
            Int32 tBinLength;
            if (tBinInput.HasValue)
            {
                if (tBinInput.Value <= 0)
                {
                    throw new ArgumentException("invalid binLength");
                }
                tBinLength = tBinInput.Value;
            }
            else
            {
                tBinLength = KTzHeader.GetPar_Int32(KTzParameters.tBin);
            }

            // preparing to run
            this.prepareToRunWithData(tBinLength);

            Console.WriteLine("" + " running for spike distribution...");

            // t = beginning of time bin, tNext = end of time bin, count = current count of spikes in consecutive time bins, countTemp = count of spikes in the current time bin
            Boolean noActivityInFullBinData;
            Int32 t = 0, count = 0, tCount = 0, countTemp = 0;
            Int32 tMax = tBinLength * (this.nSteps / tBinLength); // total time

            this.timestep = -1;
            while (t < tMax)
            {
                // running the simulation on this bin
                this.RunAndRecordTimeLength(tBinLength);

                countTemp = this.CountNonConsecutiveDataPeaks(0, tBinLength); // counts spikes only in the sampled data

                if (countTemp == 0) // found an avalanche
                {
                    if (count != 0) // the avalanche has spikes
                    {
                        this.spikeDistribution.Add(count);
                        count = 0;
                    }
                    if (tCount != 0) // and took one or more time windows
                    {
                        this.timeDistribution.Add(tCount);
                        tCount = 0;
                    }
                }
                else // keep counting the current avalanche
                {
                    count += countTemp;
                    tCount++;
                }

                noActivityInFullBinData = (countTemp == 0) && this.IsThereNoActivityInUnsampledData(); // checks the full data for spikes - needed for deciding whether to stimulate the network
                if (noActivityInFullBinData) // if there's no activity in complete sampled data
                {
                    this.ResetNeuronsAndSynapses(); // sets timestep to 0, so the stimulus will fire again
                }

                // setting the start of the next bin
                t += tBinLength;
            }
            if (count != 0) // found spikes in the last time window
            {
                this.spikeDistribution.Add(count);
            }
            if (tCount != 0) // found spikes in the last time window
            {
                this.timeDistribution.Add(tCount);
            }
        }

        /// <summary>
        /// runs this simulation in order to get the spike distribution and store it
        /// </summary>
        /// <param name="tBinInput">the time bin length to count spikes</param>
        private void RunForSpikeDistribution(Int32? tBinInput)
        {
            Int32 tBinLength;
            if (tBinInput.HasValue)
            {
                if (tBinInput.Value <= 0)
                {
                    throw new ArgumentException("invalid binLength");
                }
                tBinLength = tBinInput.Value;
            }
            else
            {
                tBinLength = KTzHeader.GetPar_Int32(KTzParameters.tBin);
            }

            // preparing to run
            this.prepareToRunWithData(tBinLength);

            Console.WriteLine("" + " running for spike distribution...");

            // t = beginning of time bin, tNext = end of time bin, count = current count of spikes in consecutive time bins, countTemp = count of spikes in the current time bin
            Boolean noActivityInCurrentData, noActivityInPreviousDataSet = false;
            Int32 t = 0, count = 0, tCount = 0, countTemp = 0;
            Int32 tMax = tBinLength * (this.nSteps / tBinLength); // total time

            if (this.samplingType == NetworkSamplingType.Full)
            {
                while (t < tMax)
                {
                    // running the simulation on this bin
                    this.timestep = t - 1;
                    this.RunAndRecordTimeLength(tBinLength);

                    countTemp = this.CountNonConsecutiveDataPeaks(0, tBinLength); // counts spikes only in the sampled data - which is the full data in this case
                    noActivityInCurrentData = countTemp == 0;//this.IsThereNoActivityInFullData(); // checks the full data for spikes - needed for deciding whether to stimulate the network

                    if (noActivityInCurrentData) // if there's no activity in current time window
                    {
                        // would waiting for two time windows make any difference in the data?
                        // I'd say that it doesn't make any difference if this.restIntervals is big
                        if (noActivityInPreviousDataSet) // and no activity in the previous time window as well
                        {
                            if (count != 0) // and if we found some spikes
                            {
                                this.spikeDistribution.Add(count); // we found an avalanche and store it
                                count = 0;
                            }
                            if (tCount != 0) // if the avalanche has one or more time windows
                            {
                                this.timeDistribution.Add(tCount); // we store avalanche duration
                                tCount = 0;
                                this.stimulus.SetNextStimulusTime(t + this.restIntervals * tBinLength); // if there's not activity for two consecutive time windows, we set the next stimulus timestep
                            }
                        }
                        else // if there's activity in the previous time window
                        {
                            noActivityInPreviousDataSet = true; // then there'll be no activity in the previous time window of the next time window
                        }
                    }
                    else // well, if there's still activity
                    {
                        count += countTemp; // we count spikes
                        tCount++; // count time windows
                        //noActivityInPreviousDataSet = false; // and, of course, there'll be activity in the previous time window of the next time window
                    }
                    // setting the start of the next bin
                    t += tBinLength;
                }
            }
            else if (this.samplingType == NetworkSamplingType.Partial)
            {
                Boolean setNextStimulus = true;
                while (t < tMax)
                {
                    // running the simulation on this bin
                    this.timestep = t - 1;
                    this.RunAndRecordTimeLength(tBinLength);

                    countTemp = this.CountNonConsecutiveDataPeaks(0, tBinLength); // counts spikes only in the sampled data
                    noActivityInCurrentData = (countTemp == 0) && this.IsThereNoActivityInUnsampledData(); // checks the full data for spikes - needed for deciding whether to stimulate the network

                    if (noActivityInCurrentData) // if there's no activity in complete sampled data
                    {
                        if (noActivityInPreviousDataSet) // and if there's no activity in previoues complete sampled data
                        {
                            if (setNextStimulus)
                            {
                                this.stimulus.SetNextStimulusTime(t + this.restIntervals * tBinLength); // set the next stimulus timestep
                                setNextStimulus = false;
                            }
                        }
                        else // if there's activity in the previous data set, but no activity in this one, then during the next time window there'll be no activity in previous data set, right?
                        {
                            noActivityInPreviousDataSet = true;
                        }
                    }
                    else
                    {
                        setNextStimulus = true;
                    }
                    //else // if there's still activity in the full data
                    //{
                    //    noActivityInPreviousDataSet = false; // then there'll be activity in the previous time window of the next time window
                    //}

                    if (countTemp == 0) // found an avalanche
                    {
                        if (count != 0) // the avalanche has spikes
                        {
                            this.spikeDistribution.Add(count);
                            count = 0;
                        }
                        if (tCount != 0) // and took one or more time windows
                        {
                            this.timeDistribution.Add(tCount);
                            tCount = 0;
                        }
                    }
                    else // keeps counting the current avalanche
                    {
                        count += countTemp;
                        tCount++;
                    }
                    // setting the start of the next bin
                    t += tBinLength;
                }
            }

            if (count != 0)
            {
                this.spikeDistribution.Add(count);
            }
            if (tCount != 0)
            {
                this.timeDistribution.Add(tCount);
            }
        }

        /// <summary>
        /// gets the output file header
        /// </summary>
        /// <returns>header for the output file</returns>
        public String GetOutputFileHeader()
        {
            String header = "";
            header += "# KT Neuron Stimulus Simulation - " + nNeurons.ToString() + " connected neurons - C#" + Environment.NewLine;
            header += "# " + DateTime.Now.ToString() + Environment.NewLine;
            header += "#-" + Environment.NewLine;
            header += KTzHeader.GetAllParamString();
            header += "#-" + Environment.NewLine;
            header += "# Initial conditions: type " + this.icType.ToString() + Environment.NewLine;
            if (icType == InitialConditionType.FromXMLFile)
            {
                header += "# " + ((GetICFromXMLFile)icGetter).fileName + Environment.NewLine;
            }
            header += "#-" + Environment.NewLine;
            header += "# Data" + Environment.NewLine;
            return header;
        }

        private String[] GetXDataColNames()
        {
            return Enumerable.Range(0, this.nNeurons).Select((i) => "x" + (i + 1).ToString()).ToArray();
        }

        private void WriteDataToFileMat(String elapsedTime)
        {
            if (hasRunned)
            {
                String oFileName = KTzHeader.GetPar_String(KTzParameters.oFile);
                if (oFileName == "") oFileName = "ktz" + this.nNeurons.ToString();

                String fileName = KTzHeader.CheckAndGetFileName(oFileName + "_sim.mat");
                String header = "# total simulation time = " + elapsedTime + Environment.NewLine + this.GetOutputFileHeader();

                var matDataBuilder = new MatFileHandler.DataBuilder();
                Dictionary<KTzParamGroups, MatFileHandler.IVariable> simParams = KTzHeader.GetAllParamPairsAsMatlabStruct(matDataBuilder);
                using (var fileStream = new System.IO.FileStream(fileName, System.IO.FileMode.Create))
                {
                    Console.WriteLine("" + " * Writing to file: {0}", fileName);
                    var writer = new MatFileHandler.MatFileWriter(fileStream);
                    writer.Write(
                        matDataBuilder.NewFile(simParams.Values.Concat(new[]{
                         matDataBuilder.NewVariable("xData", matDataBuilder.NewArray(this.xData, this.xData.Length, this.xData[0].Length)),
                         matDataBuilder.NewVariable("time", matDataBuilder.NewArray(Enumerable.Range(0, this.xData[0].Length).ToArray(), 1, this.xData[0].Length)),
                         matDataBuilder.NewVariable("spike_times", matDataBuilder.NewArray(this.spikeTimeData.ToArray(), 1, this.spikeTimeData.Count)),
                         matDataBuilder.NewVariable("spike_neurons", matDataBuilder.NewArray(this.spikeNeuronData.ToArray(), 1, this.spikeNeuronData.Count)),
                         matDataBuilder.NewVariable("file_info", matDataBuilder.NewCharArray(" xData rows -> neurons;\n xData cols -> time;\n spike_times -> time of each spike of the corresponding neuron in spike_neurons;\n spike_neurons -> neuron that spiked at the corresponding time given in spike_times")),
                         matDataBuilder.NewVariable("file_header", matDataBuilder.NewCharArray(header))
                        }))
                    );
                }

                if (this.doWriteDif)
                    this.WriteXDataDifFileTxt(oFileName + "_dif.dat", header);
                if (this.writeXDataCSVFile)
                    this.WriteDataToFileTxtInternal(oFileName + "_sim.csv", header, ",");
                // end of files
                Console.WriteLine("#-");
            }
            else
            {
                throw new Exception("The simulation hasn't runned yet!");
            }
        }

        private void WriteDataToFileTxt(String elapsedTime)
        {
            if (hasRunned)
            {
                String oFileName = KTzHeader.GetPar_String(KTzParameters.oFile);
                if (oFileName == "") oFileName = "ktz" + this.nNeurons.ToString();

                String header;
                header = "# total simulation time = " + elapsedTime + Environment.NewLine + this.GetOutputFileHeader();

                // writing dat file
                this.WriteDataToFileTxtInternal(oFileName + "_sim.dat", header, "\t");
                if (this.spikeTimeData.Count > 0)
                    Console.WriteLine("WARNING ::: spike times have been recorded but are not going to be written because output file type is Txt");
                if (this.doWriteDif)
                    this.WriteXDataDifFileTxt(oFileName + "_dif.dat", header);
                if (this.writeXDataCSVFile)
                    this.WriteDataToFileTxtInternal(oFileName + "_sim.csv", header, ",");
                // end of files
                Console.WriteLine("#-");
            }
            else
            {
                throw new Exception("The simulation hasn't runned yet!");
            }
        }

        private void WriteDataToFileTxtInternal(String fileName, String header = "", String colSep = "\t")
        {
            FileStream myFile;
            StreamWriter sw;

            fileName = KTzHeader.CheckAndGetFileName(fileName);
            Console.WriteLine("" + " * Writing to file: {0}", fileName);

            try
            {
                myFile = new FileStream(fileName, FileMode.CreateNew, FileAccess.Write);
            }
            catch (IOException)
            {
                myFile = new FileStream(fileName, FileMode.Truncate, FileAccess.Write);
            }

            sw = new StreamWriter(myFile);

            // writing header
            header += "#       t" + colSep + String.Join(colSep, this.GetXDataColNames()) + Environment.NewLine;
            sw.Write(header);

            String outputStrFmt = "{0}" + colSep + String.Join(colSep, Enumerable.Range(0, this.nNeurons).Select((n) => "{" + n.ToString() + ":0.00000000e+000}").ToArray());

            int t = 0;
            int t_total = xData[0].Length;
            while (t < t_total)
            {
                /*sw.Write("{0:0.00000000e+000}", (Double)(this.nStartStep + t));
                j = 0;
                while (j < nNeurons)
                {
                    sw.Write("\t{0:0.00000000e+000}", xData[j][t]);
                    j++;
                }
                sw.WriteLine();/**/
                sw.WriteLine(String.Format(outputStrFmt, (new object[] { (this.nStartStep + t) }).Concat(xData.GetCol(t).Cast<object>()).ToArray()));
                t++;
            }

            sw.Close();
            myFile.Close();

            Console.WriteLine("" + " * File written...");
            Console.WriteLine("" + "-");
        }

        private void WriteXDataDifFileTxt(String fileName, String header = "")
        {
            FileStream myFile;
            fileName = KTzHeader.CheckAndGetFileName(fileName);
            Console.WriteLine("" + " * Writing to file: {0}", fileName);

            try
            {
                myFile = new FileStream(fileName, FileMode.CreateNew, FileAccess.Write);
            }
            catch (IOException)
            {
                myFile = new FileStream(fileName, FileMode.Truncate, FileAccess.Write);
            }

            StreamWriter sw = new StreamWriter(myFile);
            if (header.Length > 0)
                sw.Write(header);
            sw.WriteLine("#-");
            sw.WriteLine("# Difference data between two neurons' x");
            sw.WriteLine("# t\tx(i)-x(j)");
            int m = xData[0].Length;
            int i = 0;
            int j, k;
            while (i < nNeurons)
            {
                sw.WriteLine();
                sw.WriteLine();
                sw.Write("#     t");
                j = 0;
                while (j < i)
                {
                    sw.Write("\t    x{0}-x{1}    ", i + 1, j + 1);
                    j++;
                }
                sw.WriteLine();
                k = 0;
                while (k < m)
                {
                    sw.Write("{0:0.00000000e+000}", (Double)(this.nStartStep + k));
                    j = 0;
                    while (j < i)
                    {
                        sw.Write("\t{0:0.00000000e+000}", xData[i][k] - xData[j][k]);
                        j++;
                    }
                    sw.WriteLine();
                    k++;
                }
                i++;
            }
            sw.Close();
            myFile.Close();
            Console.WriteLine("" + " * File written...");
            Console.WriteLine("" + "-");
        }


        private void GetCumulativeDist(out List<Double> normRankSize, out List<Int32> spkDistCD, out List<Double> normRankTime, out List<Int32> timeDistCD)
        {
            // we will write the distribution file as a cumulative prob. plot
            // the avalanche sizes are the first column, organized from the greatest
            // to the lowest value, while the y column is just the number of the line
            // or ocurrences greater than that value of avalanches

            // organizing the avalanche sizes
            this.spikeDistribution = this.spikeDistribution.OrderByDescending(val => val).ToList();
            this.timeDistribution = this.timeDistribution.OrderByDescending(val => val).ToList();

            //normalizing the cumulative distribution rank
            normRankSize = new List<Double>();
            spkDistCD = new List<Int32>();
            int n = this.spikeDistribution.Count;
            int i = 0;
            spkDistCD.Add(this.spikeDistribution[i]);
            normRankSize.Add((Double)(i + 1) / (Double)n);
            i++;
            while (i < n)
            {
                if (this.spikeDistribution[i] != this.spikeDistribution[i - 1])
                {
                    spkDistCD.Add(this.spikeDistribution[i]);
                    normRankSize.Add((Double)(i + 1) / (Double)n);
                }
                i++;
            }
            Double maxRank = normRankSize.Max();
            normRankSize = normRankSize.Select(d => d / maxRank).ToList();

            //normalizing the cumulative distribution rank
            normRankTime = new List<Double>();
            timeDistCD = new List<Int32>();
            n = this.timeDistribution.Count;
            i = 0;
            timeDistCD.Add(this.timeDistribution[i]);
            normRankTime.Add((Double)(i + 1) / (Double)n);
            i++;
            while (i < n)
            {
                if (this.timeDistribution[i] != this.timeDistribution[i - 1])
                {
                    timeDistCD.Add(this.timeDistribution[i]);
                    normRankTime.Add((Double)(i + 1) / (Double)n);
                }
                i++;
            }
            maxRank = normRankTime.Max();
            normRankTime = normRankTime.Select(d => d / maxRank).ToList();
        }

        private FileStream CreateOrTruncateFile(String fileName)
        {
            FileStream myFile;
            try
            {
                myFile = new FileStream(fileName, FileMode.CreateNew, FileAccess.Write);
            }
            catch (IOException)
            {
                myFile = new FileStream(fileName, FileMode.Truncate, FileAccess.Write);
            }
            return myFile;
        }

        private void WriteSpikeDistributionToFileMat(String elapsedTime, String filePrefix = "")
        {
            String oFileName = KTzHeader.GetPar_String(KTzParameters.oFile) + filePrefix;
            if (oFileName == "") oFileName = "ktz" + this.nNeurons.ToString();
            String fileName = KTzHeader.CheckAndGetFileName(oFileName + "_spkdist.mat");
            String header;

            // generating spike distributions
            if (this.spikeDistribution.Count == 0)
            {
                //throw new Exception("You should make the distribution first...");
                Console.WriteLine("WARNING: Empty distribution!");
                return;
            }

            // getting file header
            header = "# total simulation time = " + elapsedTime + Environment.NewLine + this.GetOutputFileHeader();


            List<Double> normRankSize, normRankTime;
            List<Int32> spkDistCD, timeDistCD;
            this.GetCumulativeDist(out normRankSize, out spkDistCD, out normRankTime, out timeDistCD);

            var matDataBuilder = new MatFileHandler.DataBuilder();
            Dictionary<KTzParamGroups, MatFileHandler.IVariable> simParams = KTzHeader.GetAllParamPairsAsMatlabStruct(matDataBuilder);
            using (var fileStream = new System.IO.FileStream(fileName, System.IO.FileMode.Create))
            {
                Console.WriteLine("" + " * Writing to file: {0}", fileName);
                var writer = new MatFileHandler.MatFileWriter(fileStream);
                writer.Write(
                    matDataBuilder.NewFile(simParams.Values.Concat(new[]{
                         matDataBuilder.NewVariable("s_aval", matDataBuilder.NewArray(this.spikeDistribution.ToDouble().ToArray(), 1, this.spikeDistribution.Count)),
                         matDataBuilder.NewVariable("T_aval", matDataBuilder.NewArray(this.timeDistribution.ToDouble().ToArray(), 1, this.timeDistribution.Count)),
                         matDataBuilder.NewVariable("s_CD", matDataBuilder.NewArray(spkDistCD.ToArray(), 1, spkDistCD.Count)),
                         matDataBuilder.NewVariable("T_CD", matDataBuilder.NewArray(timeDistCD.ToArray(), 1, timeDistCD.Count)),
                         matDataBuilder.NewVariable("Ps_CD", matDataBuilder.NewArray(normRankSize.ToArray(), 1, normRankSize.Count)),
                         matDataBuilder.NewVariable("PT_CD", matDataBuilder.NewArray(normRankTime.ToArray(), 1, normRankTime.Count)),
                         matDataBuilder.NewVariable("file_info", matDataBuilder.NewCharArray("s_aval -> sizes of all avalanches; T_aval -> duration of all avalanches; s_CD -> aval sizes for cumulative dist Ps_CD; Ps_CD -> cumulative dist of s_CD; T_CD -> aval duration for cumulative dist PT_CD; PT_CD -> cumulative dist of T_CD;")),
                         matDataBuilder.NewVariable("file_header", matDataBuilder.NewCharArray(header))
                    }))
                );
            }
        }

        /// <summary>
        /// get and write spike distribution to a file. We will write the distribution file as a cumulative prob plot. The avalanche sizes are the x column, organized from the greatest to the lowest value, while the y column is just the number of the line, i.e. the amount of ocurrences greater than that value of avalanches
        /// </summary>
        private void WriteSpikeDistributionToFileTxt(String elapsedTime, String filePrefix = "")
        {
            String oFileName = KTzHeader.GetPar_String(KTzParameters.oFile) + filePrefix;
            if (oFileName == "") oFileName = "ktz" + this.nNeurons.ToString();
            String header;
            FileStream myFileCD, myFileH;
            StreamWriter swCD, swH;

            // generating spike distributions
            if (this.spikeDistribution.Count == 0)
            {
                //throw new Exception("You should make the distribution first...");
                Console.WriteLine("WARNING: Empty distribution!");
                return;
            }

            // getting file header
            header = "# total simulation time = " + elapsedTime + Environment.NewLine + this.GetOutputFileHeader();


            List<Double> normRankSize, normRankTime;
            List<Int32> spkDistCD, timeDistCD;
            this.GetCumulativeDist(out normRankSize, out spkDistCD, out normRankTime, out timeDistCD);

            String sFileNameH, sFileNameCD, tFileNameH, tFileNameCD;
            GetDistTxtFileNames(oFileName, out sFileNameH, out sFileNameCD, out tFileNameH, out tFileNameCD);


            /////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////
            // Size Cumulative Distribution File
            /////////////////////////////////////////////////////////////////////////////////////

            Console.WriteLine("" + " * Writing to file: {0}", sFileNameCD);
            myFileCD = this.CreateOrTruncateFile(sFileNameCD);
            swCD = new StreamWriter(myFileCD);
            swCD.Write(header);
            swCD.WriteLine("#AvalancheSizes\tCumulativeDist");
            this.WriteTxtFileColumns(swCD, "\t", spkDistCD.ToDouble().ToArray(), normRankSize.ToArray());
            swCD.Close();
            myFileCD.Close();

            /////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////
            // Size Histogram File
            /////////////////////////////////////////////////////////////////////////////////////

            Console.WriteLine("" + " * Writing to file: {0}", sFileNameH);
            myFileH = this.CreateOrTruncateFile(sFileNameH);
            swH = new StreamWriter(myFileH);
            swH.Write(header);
            swH.WriteLine("#AvalancheSizes");
            this.WriteTxtFileColumns(swH, "\t", this.spikeDistribution.ToDouble().ToArray());
            swH.Close();
            myFileH.Close();

            /////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////
            // Time Cumulative Distribution File
            /////////////////////////////////////////////////////////////////////////////////////

            Console.WriteLine("" + " * Writing to file: {0}", tFileNameCD);
            myFileCD = this.CreateOrTruncateFile(tFileNameCD);
            swCD = new StreamWriter(myFileCD);
            swCD.Write(header);
            swCD.WriteLine("#AvalancheTimes\tCumulativeDist");
            this.WriteTxtFileColumns(swCD, "\t", timeDistCD.ToDouble().ToArray(), normRankTime.ToArray());
            swCD.Close();
            myFileCD.Close();

            /////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////
            // Time Histogram File
            /////////////////////////////////////////////////////////////////////////////////////

            Console.WriteLine("" + " * Writing to file: {0}", tFileNameH);
            myFileH = this.CreateOrTruncateFile(tFileNameH);
            swH = new StreamWriter(myFileH);
            swH.Write(header);
            swH.WriteLine("#AvalancheTimes");
            this.WriteTxtFileColumns(swH, "\t", this.timeDistribution.ToDouble().ToArray());
            swH.Close();
            myFileH.Close();


            Console.WriteLine("" + " * Files written...");
            Console.WriteLine("" + "-");
#if DEBUG
            //fileNameCD = "ktz" + nNeurons.ToString() + "n_DEBUG.csv";
            //fileNameCD = KTzHeader.CheckAndGetFileName(fileNameCD);
            //Console.WriteLine("" + " * Writing to file: {0}", fileNameCD);

            //try
            //{
            //    myFileCD = new FileStream(fileNameCD, FileMode.CreateNew, FileAccess.Write);
            //}
            //catch (IOException)
            //{
            //    myFileCD = new FileStream(fileNameCD, FileMode.Truncate, FileAccess.Write);
            //}

            //swCD = new StreamWriter(myFileCD);
            //swCD.Write("t");
            //Int32 j = 0;
            //while (j < nNeurons)
            //{
            //    swCD.Write(",x{0}", j + 1);
            //    j++;
            //}
            //swCD.WriteLine();
            //Int32 k = 0;
            //Int32 m = this.DEBUG_TData.Count;
            //while (k < m)
            //{
            //    swCD.Write("{0:0.00000000e+000}", (Double)DEBUG_TData[k]);
            //    j = 0;
            //    while (j < nNeurons)
            //    {
            //        swCD.Write(",{0:0.00000000e+000}", DEBUG_XData[j][k]);
            //        j++;
            //    }
            //    swCD.WriteLine();
            //    k++;
            //}
            //swCD.Close();
            //myFileCD.Close();
#endif
        }

        private void GetDistTxtFileNames(string oFileName, out String sFileNameH, out String sFileNameCD, out String tFileNameH, out String tFileNameCD)
        {
            sFileNameCD = oFileName + "_sCD.dat";
            sFileNameH = oFileName + "_sH.dat";
            tFileNameCD = oFileName + "_tCD.dat";
            tFileNameH = oFileName + "_tH.dat";
            sFileNameCD = KTzHeader.CheckAndGetFileName(sFileNameCD);
            sFileNameH = KTzHeader.CheckAndGetFileName(sFileNameH);
            tFileNameCD = KTzHeader.CheckAndGetFileName(tFileNameCD);
            tFileNameH = KTzHeader.CheckAndGetFileName(tFileNameH);
        }

        private void WriteTxtFileColumns(System.IO.StreamWriter sw, String colSep, params Double[][] c_data)
        {
            String outputStr = String.Join(colSep, Enumerable.Range(0, c_data.Length).Select((n) => "{" + n.ToString() + ":0.00000000e+000}").ToArray());
            int i = 0;
            //c_data = c_data.Transpose();
            while (i < c_data[0].Length)
            {
                sw.WriteLine(String.Format(outputStr, (new object[] { i }).Concat(c_data.GetCol(i).Cast<object>()).ToArray<object>()));
                i++;
            }
        }

        public Double[][] GetXData()
        {
            if (hasRunned)
            {
                return xData;
            }
            else
            {
                throw new Exception("The simulation hasn't runned yet!");
            }
        }

        public Double[] GetXData(Int32 neuronInd)
        {
            if (hasRunned)
            {
                return xData[neuronInd];
            }
            else
            {
                throw new Exception("The simulation hasn't runned yet!");
            }
        }

        //public Int32[] GetTData()
        //{
        //    if (hasRunned)
        //    {
        //        return tData;
        //    }
        //    else
        //    {
        //        throw new Exception("The simulation hasn't runned yet!");
        //    }
        //}

        /// <summary>
        /// finds the t-interval between two consecutive peaks (spikes) of the function defined by the data t,x[nIndex],
        /// assuming the peak is centered between two consecutive crossings of the x=0 axis.
        /// </summary>
        /// <param name="nIndex">the index of the neuron to find the ISI</param>
        public List<Int32> findISI(Int32 nIndex)
        {
            List<Int32> isiData = new List<Int32>();
            Int32 n = this.xData[nIndex].Length;
            Int32 i, iplus, crossCounter = 0;
            Int32 x1 = 0, x2 = 0, xs = 0, xsA = 0;

            n--;
            i = 0;
            while (this.xData[nIndex][i] >= 0.0) // this is needed to start with data that are less than zero
            {
                i++;
            }
            while (i < n)
            {
                iplus = i + 1;
                if (this.xData[nIndex][i] * this.xData[nIndex][iplus] < 0)
                {
                    if (crossCounter % 2 == 0)
                    {
                        x1 = (i + iplus) / 2;
                    }
                    else
                    {
                        x2 = (i + iplus) / 2;
                        xs = (x1 + x2) / 2;
                        isiData.Add(xs - xsA);
                        xsA = xs;
                    }
                    crossCounter++;
                }
                i++;
            }
            if (isiData.Count > 0)
            {
                isiData.RemoveAt(0); // as we start with xs = 0, the first distance calculated is now an ISI distance, but it shouldn't, so we must delete it
            }
            return isiData;
        }

        /// <summary>
        /// scans the data to get a mean bin length to divide the total time interval
        /// </summary>
        /// <param name="timeTotalLength">the total length of the time interval to divide in bins, tibically tData.Length</param>
        /// <param name="dt">the deviation from the mean</param>
        /// <returns>bin length</returns>
        private Int32 GetMeanBinLength(Int32 timeTotalLength)
        {
            List<Int32> isi = new List<Int32>();
            Int32 i = 0;
            while (i < this.nNeurons) // getting the ISI data for all the neurons
            {
                isi.AddRange(this.findISI(i));
                i++;
            }
            //Console.WriteLine(isi.Sum() / isi.Count);
            return isi.Sum() / isi.Count;// MathExtension.GreatestCommonDivisor(timeTotalLength, isi.Sum() / isi.Count);
        }

        /// <summary>
        /// divide the time in bins of defined length and counts how many spikes occur in consecutive bins, storing this value as an
        /// entry in the return list
        /// </summary>
        /// <returns>the spike ditribution</returns>
        public void MakeSpikeDistribution(Int32 tBinLength = 0)
        {
            Int32 t, n = xData[0].Length - 1, tNext, count, tCount, countTemp = 0;
            if (!this.hasRunned)
            {
                throw new Exception("simulation must run before getting spike dist!");
            }
            if (tBinLength == 0)
            {
                tBinLength = this.GetMeanBinLength(n); // 10 is the mean spike duration
            }
            else if (tBinLength < 0)// || ((n % tBinLength) != 0))
            {
                throw new ArgumentException("invalid binLength");
            }
            this.spikeDistribution = new List<Int32>();
            this.timeDistribution = new List<Int32>();
            t = 0;
            count = 0;
            tCount = 0;
            n -= tBinLength;
            while (t < n)
            {
                tNext = t + tBinLength;
                countTemp = this.CountNonConsecutiveDataPeaks(t, tNext);
                count += countTemp;
                if ((countTemp == 0) && (count != 0))
                {
                    this.spikeDistribution.Add(count);
                    this.timeDistribution.Add(tCount);
                    count = 0;
                    tCount = 0;
                }
                else
                {
                    tCount++;
                }
                t = tNext;
            }
            if (countTemp != 0)
            {
                this.spikeDistribution.Add(count);
                this.timeDistribution.Add(tCount);
            }
            //return spikeDistribution;
        }

        private void ChangeNeuronCountedState()
        {
            Int32 i, j, k, m = this.xData[0].Length, n = this.samplingIndices.Length;
            i = 0;
            while (i < n)
            {
                k = this.samplingIndices[i];
                j = 1;
                while (j < m)
                {
                    if (this.neuron[k].SpikeDetector(this.xData[k][j],this.xData[k][j-1]))//(this.xData[k][j] > 0.0D)
                    {
                        this.countedNeuron[k] = true;
                        break;
                    }
                    j++;
                }
                i++;
            }
        }

        /// <summary>
        /// counts spiking neurons only in sampled data during whole simulation time
        /// </summary>
        /// <returns>amount of neurons that have spiked during simulation</returns>
        private Int32 CountSpikingNeurons()
        {
            Int32 i, j, k, m = this.xData[0].Length, counter = 0, n = this.samplingIndices.Length;
            i = 0;
            while (i < n)
            {
                k = this.samplingIndices[i];
                j = 1;
                while (j < m)
                {
                    if (this.neuron[k].SpikeDetector(this.xData[k][j],this.xData[k][j-1]) && (!this.countedNeuron[k]))//((this.xData[k][j] > 0.0D) && (!this.countedNeuron[k]))
                    {
                        counter++;
                        break;
                    }
                    j++;
                }
                i++;
            }
            return counter;
        }

        /// <summary>
        /// check if theres any activity in the network
        /// </summary>
        /// <returns>true if theres no activity, false otherwise</returns>
        private Boolean IsThereNoActivityInUnsampledData()
        {
            if (this.unsampledIndices == null)
                return true;
            Int32 i, j, m = this.xData[this.unsampledIndices[0]].Length, n = this.unsampledIndices.Length;
            i = 0;
            while (i < n)
            {
                j = 0;
                while (j < m)
                {
                    if (this.xData[this.unsampledIndices[i]][j] > 0.0D)
                    {
                        return false;
                    }
                    j++;
                }
                i++;
            }
            return true;
        }


        /// <summary>
        /// counts all non-consecutive data peaks during all the time
        /// </summary>
        /// <param name="k">the index of the neuron</param>
        /// <returns></returns>
        public Int32 CountNonConsecutiveDataPeaks(Int32 k)
        {
            return this.CountNonConsecutiveDataPeaks(k, 0, this.xData[k].Length);
        }

        /// <summary>
        /// counts the peaks of all the neurons data between t1 and t2
        /// </summary>
        /// <param name="t1">the first index of time</param>
        /// <param name="t2">the second index of time</param>
        /// <returns>sum over the given data spikes between t1 and t2</returns>
        public Int32 CountNonConsecutiveDataPeaks(Int32 t1, Int32 t2)
        {
            Int32 peakCounter = 0;
            Int32 k = 0;//, m;
            //m = this.samplingIndices.Length;
            //k = 0;
            while (k < this.samplingN)
            {
                peakCounter += this.CountNonConsecutiveDataPeaks(this.samplingIndices[k], t1, t2);
                k++;
            }
            return peakCounter;
        }

        /// <summary>
        /// counts the peaks of the neuron k between t1 and t2
        /// </summary>
        /// <param name="k">the neuron index</param>
        /// <param name="t1">the first index of time</param>
        /// <param name="t2">the second index of time. if t2 == t1, then we count from t1 to the end of data</param>
        /// <returns>sum over k-th neuron spikes between t1 and t2</returns>
        public Int32 CountNonConsecutiveDataPeaks(Int32 k, Int32 t1, Int32 t2)
        {
            Int32 peakCounter = 0;
            t2--;
            Int32 t = t1;
            while (t < t2)
            {
                if (this.neuron[k].SpikeDetector(this.xData[k][t+1],this.xData[k][t]))
                    peakCounter++;
                t++;
            }
            return peakCounter;
        }

        /// <summary>
        /// counts the peaks of the neuron k between t1 and t2
        /// </summary>
        /// <param name="k">the neuron index</param>
        /// <param name="t1">the first index of time</param>
        /// <param name="t2">the second index of time. if t2 == t1, then we count from t1 to the end of data</param>
        /// <returns>sum over k-th neuron spikes between t1 and t2</returns>
        /*
        public Int32 CountNonConsecutiveDataPeaks(Int32 k, Int32 t1, Int32 t2)
        {
            Int32 peakCounter = 0;
            Int32 i, iNext;
            Boolean found = true;
            Boolean isNegative = (this.xData[k][t1] < 0.0);
            t2--;
            i = t1;
            while (i < t2)
            {
                iNext = i + 1;
                if (this.xData[k][i] < 0.0) // assuming that I only have a peak for y > 0
                {// if I cross the x-axis going up
                    isNegative = true;
                }
                if ((this.xData[k][i] < this.xData[k][iNext]) && (found))
                {// if I'm climbing up a peak and already found one
                    found = false; // then we should start expecting for the next peak
                }
                else if (this.xData[k][i] > 0.0)
                {// if I'm at a peak, I only count it if it is greater than 0
                    if ((this.xData[k][i] > this.xData[k][iNext]) && (!found) && (isNegative))
                    {// I'm only at a peak if the next point is smaller than the current one, if I didn't already found a peak and if I have previously crossed the x-axis going up
                        found = true;
                        isNegative = false;
                        peakCounter++;
                    }
                }
                i++;
            }
            return peakCounter;
        }/**/

    }
}
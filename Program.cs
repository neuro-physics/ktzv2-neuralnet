using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using KTzV2.Data;
using KTzV2.Data.Header;
using KTzV2.Sims.Network;
using MatFileHandler;

namespace KTzV2
{
    public enum SimulationType
    {
        Dynamics,
        DynamicsWithinParamRange,
        Bifurcation,
        KTzPhaseDiagram
    }

    public enum ParamForRangeInDynamicsSim
    {
        J, I, alpha, r, Theta
    }

    public enum DynamicsSimType
    {
        ContinuousTime,
        NetworkReset
    }

    public enum BifurcationWritePolicy
    {
        OnTheFly,
        OnTheEnd
    }

    public enum OutputAverageMode
    {
        OverRealizations, // takes average over n independent avalanches
        OverTime // takes average over temporal fluctuations of activated density
    }

    class Program
    {
        /*static void Main(string[] args)
        {
            double[] y = new double[1000];
            double dx = 0.01;
            for (int i = 0; i < 1000; i++)
            {
                y[i] = Math.Cos((Double)i * dx);
            }

            double[] dy1 = (new KTzV2.Maths.Calculus.FirstOrderDerivative(dx, y, KTzV2.Maths.Calculus.DerivativeMethod.CentralDiff)).GetDerivative();
            double[] dy2 = (new KTzV2.Maths.Calculus.FirstOrderDerivative(dx, y, KTzV2.Maths.Calculus.DerivativeMethod.ImprovedCentralDiff)).GetDerivative();
            double[] dy3 = (new KTzV2.Maths.Calculus.FirstOrderDerivative(dx, y, KTzV2.Maths.Calculus.DerivativeMethod.RichardsonExtrap)).GetDerivative();

            System.IO.FileStream fs = new System.IO.FileStream("teste.dat", System.IO.FileMode.Create, System.IO.FileAccess.Write);
            System.IO.StreamWriter sw = new System.IO.StreamWriter(fs);
            sw.WriteLine("#x\ty\tdy_Cent\tdy_ImpC\tdy_Rich");
            for (int i = 0; i < 1000; i++)
            {
                sw.WriteLine("{0:0.00000000e+000}\t{1:0.00000000e+000}\t{2:0.00000000e+000}\t{3:0.00000000e+000}\t{4:0.00000000e+000}", (Double)i * dx, y[i], dy1[i], dy2[i], dy3[i]);
            }
            sw.Close();
            fs.Close();
        }*/

        /// <summary>
        /// Initiates data variables for Bifurcation simulation. Returns variable[nPar2][nPar1]
        /// </summary>
        /// <param name="nPar2">amount of values in first array dimension</param>
        /// <param name="nPar1">amount of values in second array dimension</param>
        static void InitiateBifurcationDataVariables(int nPar2, int nPar1, out Double[][] mMax, out Double[][] m, out Double[][] mStdDev, out Double[][] chi, out Double[][] chiStdDev, out Double[][] U4, out Double[][] U4StdDev, out Double[][] dU4, out Double[][] dU4Err)
        {
            mMax = new Double[nPar2][];
            m = new Double[nPar2][];
            mStdDev = new Double[nPar2][];
            chi = new Double[nPar2][];
            chiStdDev = new Double[nPar2][];
            U4 = new Double[nPar2][];
            U4StdDev = new Double[nPar2][];
            dU4 = new Double[nPar2][];
            dU4Err = new Double[nPar2][];

            // the increment of the I and J variables
            //Double dPar1 = (nPar1 > 1 ? (maxPar1 - minPar1) / Convert.ToDouble(nPar1 - 1) : 0.0);
            //Double dPar2 = (nPar2 > 1 ? (maxPar2 - minPar2) / Convert.ToDouble(nPar2 - 1) : 0.0);

            Int32 i, j;

            // resetting the values of all the data variables
            i = 0;
            while (i < m.Length)
            {
                mMax[i] = new Double[nPar1];
                m[i] = new Double[nPar1];
                mStdDev[i] = new Double[nPar1];
                chi[i] = new Double[nPar1];
                chiStdDev[i] = new Double[nPar1];
                U4[i] = new Double[nPar1];
                U4StdDev[i] = new Double[nPar1];
                dU4[i] = new Double[nPar1];
                dU4Err[i] = new Double[nPar1];

                j = 0;
                while (j < m[i].Length)
                {
                    mMax[i][j] = 0.0D;
                    m[i][j] = 0.0D;
                    mStdDev[i][j] = 0.0D;
                    chi[i][j] = 0.0D;
                    chiStdDev[i][j] = 0.0D;
                    U4[i][j] = 0.0D;
                    U4StdDev[i][j] = 0.0D;
                    dU4[i][j] = 0.0D;
                    dU4Err[i][j] = 0.0D;

                    j++;
                }
                i++;
            }
        }

        static void Main(string[] args)
        {
#if DEBUG
            /**
            var colNames = new String[] { "d1", "d2" };
            var d1 = new Double[] { 1, 2, 3 };
            var d2 = new Double[] { 4, 5, 6 };
            System.IO.StreamWriter mysw = KTzV2.Data.Header.KTzHeader.CreateOutTxtFile("test.txt", "this is my header", new string[] { });
            KTzV2.Data.Header.KTzHeader.WriteTxtFileColumns(mysw, colNames, d1, d2);
            KTzV2.Data.Header.KTzHeader.CloseOutTxtFile(mysw);
            return;
            /**/
            /*
             * testing basic matlab writing capabilities
             * /
            var b = Enumerable.Range(0,10000).Select(el => Convert.ToDouble(el)).ToArray();
            var a = new Double[][] { b, b };
            var matDataBuilder = new MatFileHandler.DataBuilder();
            var array = matDataBuilder.NewArray<Double>(a, a.Length, a[0].Length);
            var variable = matDataBuilder.NewVariable("a", array);
            var matFile = matDataBuilder.NewFile(new[] { variable });
            using (var fileStream = new System.IO.FileStream("test.mat", System.IO.FileMode.Create))
            {
                var writer = new MatFileHandler.MatFileWriter(fileStream);
                writer.Write(matFile);
                Console.WriteLine("done");
            }
            return;
            /**/
            /* testing advanced matlab writing capabilities
             *
            KTzV2.Data.Header.KTzHeader.ResetAllParamLists();
            var p1 = new Double[] { 1, 2, 3 };
            var p2 = new Double[] { 4, 5, 6, 7 };
            var z0 = new Double[][] { new Double[] { 1, 2, 3, 4 }, new Double[] { 5, 6, 7, 8 }, new Double[] { 9, 10, 11, 12 } };
            KTzV2.Data.Header.KTzHeader.Write11ColFileMat(ref p1, ref p2,
            ref z0, ref z0, ref z0, ref z0, ref z0, ref z0, ref z0, ref z0, ref z0,
            new String[] {"c1","c2","c3","c4","c5","c6","c7","c8","c9","c10","c11" }, "test", "my file header");
            Console.WriteLine("Press key to continue...");
            Console.ReadKey();
            return;
            /**/
#endif
            KTzHeader.ResetAllParamLists();
            if (!KTzHeader.CheckForInputArguments(args)) return;
#if !DEBUG
            try
            {
#endif
                System.Diagnostics.Stopwatch chronometer = new System.Diagnostics.Stopwatch();
                chronometer.Start();

                if (((SimulationType)KTzHeader.GetPar_Int32(KTzParameters.simType) == SimulationType.Dynamics) || ((SimulationType)KTzHeader.GetPar_Int32(KTzParameters.simType) == SimulationType.DynamicsWithinParamRange) )
                {
                    KTzNetworkSimulator KTzProg = new KTzNetworkSimulator(false);

                    var parForRange = (ParamForRangeInDynamicsSim)KTzHeader.GetPar_Int32(KTzParameters.ParamForRange);
                    bool hasRangeParam = Enum.TryParse(parForRange.ToString() + "Range", false, out KTzParameters _par);
                    if (!hasRangeParam)
                        throw new ArgumentException("Attempted to create a Range for a parameter that has no Range input: " + parForRange.ToString());


                    //KTzParameters par = KTzParameters.J;
                    Enum.TryParse(parForRange.ToString(), false, out KTzParameters par);
                    String parName = parForRange.ToString();
                    Double[] parRange;

                    if ((SimulationType)KTzHeader.GetPar_Int32(KTzParameters.simType) == SimulationType.DynamicsWithinParamRange)
                        parRange = KTzHeader.GetRangeFor(par);
                    else
                        parRange = new double[] { KTzHeader.GetPar_Double(par) };


                    KTzV2.Neurons.NeuronType nt = (KTzV2.Neurons.NeuronType)KTzHeader.GetPar_Int32(KTzParameters.neuron);
                    if ((nt == KTzV2.Neurons.NeuronType.SIElement) && (parForRange != ParamForRangeInDynamicsSim.Theta))
                        Console.WriteLine("WARNING: SIElement is the neuron, but ParamForRange is not Theta");
                    //{
                        //par = KTzParameters.Theta;
                        //parRange = KTzHeader.GetRangeFor(KTzParameters.Theta);
                        //parName = "Theta";
                    //}

                    if (parRange.Length > 1)
                    {
                        if (nt == KTzV2.Neurons.NeuronType.SIElement)
                            Console.WriteLine("WARNING: running Avalanche/Dynamics simulation for many Theta defined by [minTheta,maxTheta] or ThetaRange!");
                        else
                            Console.WriteLine("WARNING: running Avalanche/Dynamics simulation for many J defined by [minJ,maxJ] or JRange!");
                    }

                    if (KTzHeader.GetPar_Int32(KTzParameters.wSpk) == 1)
                    {
                        chronometer.Reset();
                        Int32 k = 0;
                        foreach (Double parValue in parRange)
                        {
                            KTzProg.SetParam(par, parValue);
                            Console.WriteLine("*** (sim {0}/{1}) Running for Avalanches {3} = {2}", ++k, parRange.Length, parValue, parName);
                            chronometer.Start();
                            KTzProg.RunForAvalanches(KTzHeader.GetPar_Int32(KTzParameters.tBin));
                            chronometer.Stop();
                            KTzProg.WriteSpikeDistributionToFile(chronometer.Elapsed.ToString(), String.Format("_{1}{0}", parValue, parName));
                            chronometer.Reset();
                        }
                    }
                    if (KTzHeader.GetPar_Int32(KTzParameters.wData) == 1)
                    {
                        chronometer.Reset();
                        Int32 k = 0;
                        foreach (Double parValue in parRange)
                        {
                            KTzProg.SetParam(par, parValue);
                            Console.WriteLine("*** (sim {0}/{1}) Running for Dynamics {3} = {2}", ++k, parRange.Length, parValue, parName);
                            chronometer.Start();
                            Int32? iStim = KTzHeader.GetPar_Int32(KTzParameters.iStim);
                            if (iStim.Value == -1)
                                iStim = null;
                            KTzProg.Run(iStim);
                            chronometer.Stop();
                            KTzProg.WriteDataToFile(chronometer.Elapsed.ToString());
                            chronometer.Reset();
                        }
                    }
                }
                else if ((SimulationType)KTzHeader.GetPar_Int32(KTzParameters.simType) == SimulationType.Bifurcation)
                {
                    BifurcationWritePolicy wPol = (BifurcationWritePolicy)KTzHeader.GetPar_Int32(KTzParameters.bifWrite);

                    String filePrefix = KTzHeader.GetPar_String(KTzParameters.oFile);
                    String oFileName = filePrefix;
                    if (filePrefix == "") oFileName = "sb";

                    if (KTzHeader.GetPar_Int32(KTzParameters.outAvgMode) == (Int32)KTzV2.OutputAverageMode.OverTime)
                    {
                        KTzHeader.SetPar(KTzParameters.dynType, (Int32)KTzV2.DynamicsSimType.ContinuousTime);
                    }

                    //Int32 nNeurons = Convert.ToInt32(Math.Pow(KTzHeader.GetPar_Int32(KTzParameters.Lx), KTzHeader.GetPar_Int32(KTzParameters.dim)));
                    KTzNetworkSimulator KTzProg = new KTzNetworkSimulator(false);
                    Int32 nNeurons = KTzProg.nNeurons;

                    Boolean wAvalObs = KTzHeader.GetPar_Int32(KTzParameters.wObs) == 1;

                    KTzParameters par1, par2;

                    // selecting the primary parameter
                    KTzV2.Synapses.SynapseType st = (KTzV2.Synapses.SynapseType)KTzHeader.GetPar_Int32(KTzParameters.sType);
                    if (st == KTzV2.Synapses.SynapseType.KTDynamicChemicalSynapse)
                    {
                        Console.WriteLine("PRIMARY BIFURCATION PARAM: alpha since sType is KTDynamicChemicalSynapse");
                        par1 = KTzParameters.alpha;
                    }
                    else
                    {
                        Console.WriteLine("PRIMARY BIFURCATION PARAM: J since sType is not KTDynamicChemicalSynapse");
                        par1 = KTzParameters.J;
                    }

                    KTzV2.Neurons.NeuronType nt = (KTzV2.Neurons.NeuronType)KTzHeader.GetPar_Int32(KTzParameters.neuron);
                    if (nt == KTzV2.Neurons.NeuronType.SIElement)
                    {
                        Console.WriteLine("WARNING: forcing PRIMARY BIFURCATION PARAM == Theta since neuron is SIElement");
                        par1 = KTzParameters.Theta;
                    }

                    // selecting the secondary parameter
                    KTzV2.Stimuli.StimulusType stt = (KTzV2.Stimuli.StimulusType)KTzHeader.GetPar_Int32(KTzParameters.stimType);
                    if (stt == KTzV2.Stimuli.StimulusType.PoissonProcess)
                    {
                        Console.WriteLine("SECONDARY BIFURCATION PARAM: r since stimType is PoissonProcess");
                        par2 = KTzParameters.r;
                    }
                    else
                    {
                        Console.WriteLine("SECONDARY BIFURCATION PARAM: I since stimType is not PoissonProcess");
                        par2 = KTzParameters.I;
                    }
                    
                    // getting the setted parameters
                    String par1Name = par1.ToString();
                    String par2Name = par2.ToString();
                    
                    // selecting variable to count
                    KTzV2.Data.CountVariable cVar = (KTzV2.Data.CountVariable)KTzHeader.GetPar_Int32(KTzParameters.cVar);

                    Int32 timeWindow = KTzHeader.GetPar_Int32(KTzParameters.tBin);

                    Int32 nSim = KTzHeader.GetPar_Int32(KTzParameters.nSim);

                    if ((stt == KTzV2.Stimuli.StimulusType.Delta) || (stt == KTzV2.Stimuli.StimulusType.DeltaWhenInactive))
                    {
                        Console.WriteLine("WARNING: the program will run {0} independent avalanches since stimType is Delta or DeltaWhenInactive", nSim);
                    }
                    else
                    {
                        if (cVar == KTzV2.Data.CountVariable.NumberOfNeurons)
                        {
                            Console.WriteLine("WARNING: forcing count of spikes since stimType is not Delta nor DeltaWhenInactive");
                            cVar = KTzV2.Data.CountVariable.NumberOfSpikes;
                            KTzHeader.SetPar(KTzParameters.cVar, (int)cVar);
                        }
                        Console.WriteLine("WARNING: each of the {0} simulations will run for {1} time steps since stimType is not Delta nor DeltaWhenInactive", nSim, KTzHeader.GetPar_Int32(KTzParameters.nSteps));
                    }

                    if ((stt == KTzV2.Stimuli.StimulusType.PoissonProcess) && (cVar == KTzV2.Data.CountVariable.NumberOfNeurons))
                    {
                        Console.WriteLine("WARNING: forcing count of spikes since stimType is PoissonProcess");
                        cVar = KTzV2.Data.CountVariable.NumberOfSpikes;
                        KTzHeader.SetPar(KTzParameters.cVar, (int)cVar);
                    }

                    // set of J and I to use in the simulations
                    Double[] Par1 = KTzHeader.GetRangeFor(par1);//new Double[nPar1]; // +1 to include the Jmax as well
                    Double[] Par2 = KTzHeader.GetRangeFor(par2);//new Double[nPar2]; // +1 to include the Imax as well

                    Int32 nPar1 = Par1.Length;
                    Int32 nPar2 = Par2.Length;

                    // data variables
                    Int32[] countPerSim = new Int32[nSim]; // total number of spikes in one simulation used to calculate the standard deviation
                    ThermoStatistics ts;

                    // the increment of the I and J variables
                    //Double dPar1 = (nPar1 > 1 ? (maxPar1 - minPar1) / Convert.ToDouble(nPar1 - 1) : 0.0);
                    //Double dPar2 = (nPar2 > 1 ? (maxPar2 - minPar2) / Convert.ToDouble(nPar2 - 1) : 0.0);
                    
                    Int32 i, j;

                    // resetting the values of all the data variables
                    Double[][] mMax = new Double[nPar2][];
                    Double[][] m = new Double[nPar2][];
                    Double[][] mStdDev = new Double[nPar2][];
                    Double[][] chi = new Double[nPar2][];
                    Double[][] chiStdDev = new Double[nPar2][];
                    Double[][] U4 = new Double[nPar2][];
                    Double[][] U4StdDev = new Double[nPar2][];
                    Double[][] dU4 = new Double[nPar2][];
                    Double[][] dU4Err = new Double[nPar2][];/**/
                    InitiateBifurcationDataVariables(nPar2, nPar1, out mMax, out m, out mStdDev, out chi, out chiStdDev, out U4, out U4StdDev, out dU4, out dU4Err);

                    // adjusting the output file header
                    String fHeader;
                    fHeader =  "# Network response for stimulus by varying coupling parameter and stimulus intensity (a try of phase (bifurcation) diagram)" + Environment.NewLine;
                    fHeader += "#-" + Environment.NewLine;
                    fHeader += KTzProg.GetOutputFileHeader();
                    fHeader += "#-" + Environment.NewLine;
                    fHeader += "# " + par1 + " = as expressed in the data below" + Environment.NewLine;
                    fHeader += "# " + par2 + " = as expressed in the data below" + Environment.NewLine;
                    fHeader += "# total simulation steps = " + KTzHeader.GetPar_Int32(KTzParameters.nSteps).ToString();
                    String varName = cVar == CountVariable.NumberOfNeurons ? "nNeu" : "nSpk";
                    String[] oFileColNames = new String[] { par2Name, par1Name, varName + "_max", varName + "_M", varName + "_SD", "chi_M", "chi_SD", "U4_M", "U4_SD", "dU4", "dU4_Err" };
                    String outputStrFmt = String.Join("\t", Enumerable.Range(0, oFileColNames.Length).Select((x) => "{" + Convert.ToString(x) + ":0.00000000e+000}").ToArray());

                    System.IO.StreamWriter sw;

                    String tmpFileName = oFileName + "_tmp.dat";
                    if (wPol == BifurcationWritePolicy.OnTheFly)
                    {
                        sw = KTzHeader.CreateOutTxtFile(tmpFileName, fHeader, oFileColNames);
                    }
                    else
                    {
                        sw = null;
                    }

                    chronometer.Reset();
                    chronometer.Start();
                    Int32 simCounter = 0, simTotal = nPar2 * nPar1;
                    i = 0;
                    while (i < Par2.Length) // stimulus loop
                    {
                        // calculating the I of the current simulation
                        //Par2[i] = minPar2 + (Double)i * dPar2;
                        KTzProg.SetParam(par2, Par2[i]);

                        j = 0;
                        while (j < Par1.Length) // synaptic coupling loop
                        {
                            simCounter++;

                            // calculating the J of the current simulation
                            //Par1[j] = minPar1 + (Double)j * dPar1;
                            KTzProg.SetParam(par1, Par1[j]);
                            
                            Console.WriteLine("*** (sim {4}/{5}) Running ({0},{1}) = ({2},{3})", par1Name, par2Name, Par1[j], Par2[i], simCounter, simTotal);

                            if (KTzHeader.GetPar_Int32(KTzParameters.outAvgMode) == (Int32)KTzV2.OutputAverageMode.OverRealizations)
                            {
                                if ((stt == KTzV2.Stimuli.StimulusType.Delta) || (stt == KTzV2.Stimuli.StimulusType.DeltaWhenInactive))
                                {
                                    ts = KTzProg.RunNAvalanches();
                                }
                                else
                                {
                                    ts = KTzProg.RunNTimesForSpikeCount();
                                }
                            }
                            else //(KTzHeader.GetPar_Int32(KTzParameters.outAvgMode) == (Int32)KTz.OutputAverageMode.OverTime)
                            {
                                ts = KTzProg.RunForSpikeCountTemporalDynamics();
                            }
                            if (wAvalObs)
                            {
                                ts.WriteObservations(KTzHeader.CheckAndGetFileName(String.Format("avs_{4}_{0}{1}_{2}{3}", par1Name, Par1[j], par2Name, Par2[i], filePrefix)),
                                    "# avalanches for the following simulation" + Environment.NewLine + KTzProg.GetOutputFileHeader() + "#aval_size", (KTzV2.Data.Header.OutputFileFormat)KTzV2.Data.Header.KTzHeader.GetPar_Int32(KTzParameters.oFileFormat));
                            }
                            mMax[i][j] = ts.mMax;
                            m[i][j] = ts.m;
                            mStdDev[i][j] = ts.mStdDev;
                            chi[i][j] = ts.chi;
                            chiStdDev[i][j] = ts.chiStdDev;
                            U4[i][j] = ts.U4;
                            U4StdDev[i][j] = ts.U4StdDev;

                            if (wPol == BifurcationWritePolicy.OnTheFly)
                            {
                                // String[] oFileColNames = new String[] { par2Name, par1Name, varName + "_max", varName + "_M", varName + "_SD", "chi_M", "chi_SD", "U4_M", "U4_SD", "dU4", "dU4_Err" };
                                sw.WriteLine(outputStrFmt, Par2[i], Par1[j], mMax[i][j], m[i][j], mStdDev[i][j], chi[i][j], chiStdDev[i][j], U4[i][j], U4StdDev[i][j], Double.NaN, Double.NaN);
                                // write to opened file here
                                // dU4 == NaN
                            }

                            // going to a new simulation
                            j++;
                        }
                        dU4[i] = (new KTzV2.Maths.Calculus.FirstOrderDerivative(Par1, U4[i], KTzV2.Maths.Calculus.DerivativeMethod.ImprovedCentralDiff)).GetDerivative(out dU4Err[i]);
                        i++;
                    }

                    if (wPol == BifurcationWritePolicy.OnTheFly)
                    {
                        KTzHeader.CloseOutTxtFile(sw);
                    }

                    chronometer.Stop(); // measuring simulation total time
                    fHeader += "#-" + Environment.NewLine;
                    fHeader += "# total simulation time: " + chronometer.Elapsed.ToString() + Environment.NewLine + "#-";

                    KTzHeader.Write11ColFile<Double, Double>(ref Par1, ref Par2,
                            ref mMax, ref m, ref mStdDev, ref chi, ref chiStdDev, ref U4, ref U4StdDev, ref dU4, ref dU4Err,
                            oFileColNames,
                            oFileName, fHeader);
                    if ((System.IO.File.Exists(tmpFileName)) && (wPol == BifurcationWritePolicy.OnTheFly))
                    {
                        System.IO.File.Delete(tmpFileName);
                    }
                }
                else if ((SimulationType)KTzHeader.GetPar_Int32(KTzParameters.simType) == SimulationType.KTzPhaseDiagram)
                {
                    Console.WriteLine("Nothing to do here!");
                    throw new NotImplementedException();
                    // calculates the phase diagram for one single KTz neuron
                }
           
#if !DEBUG
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
            }
#endif
            if (KTzHeader.WAIT_ON_FINISH)
            {
                Console.WriteLine("Finished! Press any key to exit...");
                //Console.Beep(2000, 500); // beeeeeep
                Console.ReadKey();
            }
        }
    }
}

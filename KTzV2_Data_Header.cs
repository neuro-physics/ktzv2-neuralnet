using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MatFileHandler;
/*
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.Data.Matlab;
/**/
namespace KTzV2.Data.Header
{
    public enum KTzParamGroups
    {
        Network,
        Neuron,
        InitCond,
        Synapse,
        Stimulus,
        Simulation
    }

    public enum OutputFileFormat
    {
        txt,
        mat
    }

    public enum QuenchedDisorderType
    {
        Gaussian,
        Uniform,
        None
    }

    public class QuenchedDisorderParse
    {
        public KTzParameters Param { get; private set; }
        public QuenchedDisorderType Type { get; private set; }
        public Double Mean { get; private set; }
        public Double Stddev { get; private set; }
        public Double Min { get; private set; }
        public Double Max { get; private set; }
        private KTzV2.Maths.Random.GaussianRand GRand { get; set; }
        private KTzV2.Maths.Random.HomogeneousRand URand { get; set; }
        public Func<Double> NextSample { get; private set; }
        public Action<Double> Set { get; private set; }
        public QuenchedDisorderParse(KTzParameters param, String setting)
        {
            this.Param = param;
            this.Parse(setting);
        }
        private void Parse(String setting)
        {
            String[] ss = setting.Split(':');
            if (ss.Length != 3)
                throw new ArgumentException(String.Format("Disorder setting for {0} must have two numbers: {0}Disorder=Uniform_or_Gaussian:min_or_mean:max_or_stddev",this.Param));
            try
            {
                this.Type = (QuenchedDisorderType)Enum.Parse(typeof(QuenchedDisorderType), ss[0], true);
            }
            catch (ArgumentException)
            {
                Console.WriteLine("WARNING: Unknown disorder type for {0}: {1}", this.Param, ss[0]);
            }
            switch (this.Type)
            {
                case QuenchedDisorderType.Gaussian:
                    this.Mean = Convert.ToDouble(ss[1]);
                    this.Stddev = Convert.ToDouble(ss[2]);
                    this.Min = Double.NegativeInfinity;
                    this.Max = Double.PositiveInfinity;
                    this.GRand = new KTzV2.Maths.Random.GaussianRand(this.Mean, this.Stddev);
                    this.NextSample = this.NextSample_G;
                    this.Set = this.SetG;
                    break;
                case QuenchedDisorderType.Uniform:
                    this.Min = Convert.ToDouble(ss[1]);
                    this.Max = Convert.ToDouble(ss[2]);
                    this.Mean = (this.Min + this.Max) / 2.0D;
                    this.Stddev = Math.Sqrt( Math.Pow(this.Max - this.Min,2.0D)/12.0D );
                    this.URand = new KTzV2.Maths.Random.HomogeneousRand();
                    this.NextSample = this.NextSample_U;
                    this.Set = this.SetU;
                    break;
                case QuenchedDisorderType.None:
                    this.SetNone(Convert.ToDouble(ss[1]));
                    this.NextSample = () => this.Mean;
                    this.Set = this.SetNone;
                    break;
                default:
                    throw new ArgumentException("unknown disorder type");
            }
        }

        private void SetNone(Double Mean)
        {
            this.Mean = Mean;
            this.Stddev = 0.0D;
            this.Min = this.Mean;
            this.Max = this.Mean;
        }

        private void SetU(Double Mean)
        {
            var d = Math.Abs(this.Max - this.Min);
            this.Mean = Mean;
            this.Min = this.Mean - d/2.0D;
            this.Max = this.Mean + d/2.0D;
            this.Stddev = Math.Sqrt(Math.Pow(this.Max - this.Min, 2.0D) / 12.0D);
        }

        private void SetG(Double Mean)
        {
            Double cv = this.Stddev / this.Mean;
            this.Mean = Mean;
            this.Stddev = cv * this.Mean;
            this.GRand = new KTzV2.Maths.Random.GaussianRand(this.Mean, this.Stddev);
        }

        public override String ToString()
        {
            if (this.Type != QuenchedDisorderType.None)
                return String.Format("{0}_Disorder={1}:(min,max,mean,std)=({2},{3},{4},{5})",this.Param,this.Type,this.Min,this.Max,this.Mean,this.Stddev);
            return String.Format("{0}_Disorder={1}",this.Param,this.Type);
        }

        private Double NextSample_U()
        {
            return (this.Max - this.Min) * this.URand.GetRandomFull() + this.Min;
        }

        private Double NextSample_G()
        {
            return this.GRand.getRandom();
        }
    }

    public static class MyExtensionMethods
    {
        public static MatFileHandler.IArrayOf<T> NewArray<T>(this MatFileHandler.DataBuilder matDataBuilder, T[][] data, params int[] dimensions)
            where T : struct
        {
            var array = matDataBuilder.NewArray<T>(dimensions);
            var m = data.Length;
            var n = data[0].Length;
            for (var i = 0; i < m; i++)
            {
                for (var j = 0; j < n; j++)
                {
                    array[j * m + i] = data[i][j];
                    //array[i * n + j] = data[i][j];
                }
            }
            return array;
        }

        public static T[][] Transpose<T>(this T[][] arr)
            where T : struct
        {
            int rowCount = arr.Length;
            int columnCount = arr[0].Length;
            T[][] transposed = new T[columnCount][];
            if (rowCount == columnCount)
            {
                transposed = (T[][])arr.Clone();
                for (int i = 1; i < rowCount; i++)
                {
                    for (int j = 0; j < i; j++)
                    {
                        T temp = transposed[i][j];
                        transposed[i][j] = transposed[j][i];
                        transposed[j][i] = temp;
                    }
                }
            }
            else
            {
                for (int column = 0; column < columnCount; column++)
                {
                    transposed[column] = new T[rowCount];
                    for (int row = 0; row < rowCount; row++)
                    {
                        transposed[column][row] = arr[row][column];
                    }
                }
            }
            return transposed;
        }

        public static T[] GetCol<T>(this T[][] X, int j)
            where T: struct
        {
            T[] c = new T[X.Length];
            int i = 0;
            while (i < X.Length)
            {
                c[i] = X[i][j];
                i++;
            }
            return c;
        }

        public static IEnumerable<Double> ToDouble(this List<int> x)
        {
            return x.Select(el => Convert.ToDouble(el));
        }
    }

    public static class KTzHeader
    {
        public static Boolean WAIT_ON_FINISH = false;
        public static Dictionary<KTzParameters, Double> ParamList_Double { get; private set; }
        public static Dictionary<KTzParameters, Int32> ParamList_Int32 { get; private set; }
        public static Dictionary<KTzParameters, String> ParamList_String { get; private set; }
        public static Dictionary<KTzParameters, String> ParamDescription { get; private set; }
        public static Dictionary<KTzParameters, KTzParamGroups> ParamGroups { get; private set; }
        public static Dictionary<KTzParameters, Func<Int32, String>> ParamConverter { get; private set; }
        public static HashSet<KTzParameters> ParamAllowedQuenchedDisorder { get; private set; }
        public static Dictionary<KTzParameters, QuenchedDisorderParse> ParamDisorderSetting { get; private set; }
        public static String inputArgsStr { get; private set; }

        public static bool CanHaveDisorder(KTzParameters p, bool silent = true)
        {
            bool hasDisorderParam = Enum.TryParse(p.ToString() + "Disorder", false, out KTzParameters pdis);
            if ((!hasDisorderParam) && (!silent))
                Console.WriteLine("WARNING: {0}Disorder is not defined in the KTzParameters enum and its disorder cannot be set", p);
            return hasDisorderParam && ParamAllowedQuenchedDisorder.Contains(p);
        }
        public static bool CanHaveDisorder(String p, bool silent = true)
        {
            try
            {
                return KTzHeader.CanHaveDisorder((KTzParameters)Enum.Parse(typeof(KTzParameters), p, false), silent: silent);
            }
            catch (ArgumentException)
            {
                Console.WriteLine("WARNING: {0} is not a valid param name", p);
                return false;
            }
        }
        public static bool IsDisorderInputParam(KTzParameters p, bool silent = true)
        {
            return KTzHeader.IsDisorderInputParam(p.ToString(), silent: silent);
        }
        public static bool IsDisorderInputParam(String p, bool silent=true)
        {
            String par = p.Replace("Disorder", "");
            if ((!KTzHeader.CanHaveDisorder(par)) && (!silent))
                Console.WriteLine("WARNING: {0} cannot have quenched disorder", par);
            return p.Contains("Disorder") && KTzHeader.CanHaveDisorder(par, silent: silent);
        }

        public static QuenchedDisorderParse GetDisorderSetting(KTzParameters p)
        {
            return KTzHeader.ParamDisorderSetting[p];
        }

        public static Boolean CheckForInputArguments(String[] args)
        {
            if ((args.Length == 0) || (!args.Contains<String>("-run")))
            {
                KTzHeader.ShowHelp();
                return false;
            }

            String temp = args[0].ToLower();
            if (temp[0] == '-') temp = temp.Substring(1);
            if (temp == "help")
            {
                KTzHeader.ShowHelp();
                return false;
            }

            if (args.Contains("-wait")) KTzHeader.WAIT_ON_FINISH = true;
            KTzParameters p;
            KTzHeader.inputArgsStr = KTzHeader.GetKTzV2ExeName();
            foreach (String str in args)
            {
                KTzHeader.inputArgsStr += " " + str;
                if ((str == "-ask") || (str == "-run") || (str == "-wait")) continue;
                String[] parVal = str.Split(new string[] { "=" }, StringSplitOptions.RemoveEmptyEntries);
                if (parVal[0][0] == '-') parVal[0] = parVal[0].Substring(1); // ignoring the "-" in the beginning
                try
                {
                    p = (KTzParameters)Enum.Parse(typeof(KTzParameters), parVal[0], false);
                }
                catch (ArgumentException)
                {
                    KTzHeader.ErrorUnknownParameter(parVal[0], parVal[1]);
                    return false;
                }
                try
                {
                    if (KTzHeader.ParamList_Double.ContainsKey(p))
                    {
                        KTzHeader.SetPar(p, Convert.ToDouble(parVal[1]));
                    }
                    else
                    {
                        KTzHeader.SetPar(p, Convert.ToInt32(parVal[1]));
                    }
                }
                catch (FormatException)
                {
                    try
                    {
                        KTzHeader.SetParByText(p, parVal[1]);
                    }
                    catch (ArgumentOutOfRangeException e)
                    {
                        KTzHeader.ErrorUnacceptedValue(parVal[0], parVal[1]);
                        Console.WriteLine(e.Message);
                        return false;
                    }
                    catch (ArgumentException e)
                    {
                        KTzHeader.ErrorUnknownParameter(parVal[0], parVal[1]);
                        Console.WriteLine(e.Message);
                        return false;
                    }
                }
                catch (ArgumentException e)
                {
                    KTzHeader.ErrorUnknownParameter(parVal[0], parVal[1]);
                    Console.WriteLine(e.Message);
                    return false;
                }
                Console.WriteLine("Set: {0} = {1}", parVal[0], parVal[1]);
            }

            if (args.Contains<String>("-ask"))
            {
                Console.Write(KTzHeader.GetAllParamString(true));
                KTzHeader.AskForParamsChanges();
            }

            return true;
        }

        private static void SetParByText(KTzParameters p, String value)
        {
            Object val;
            switch (p)
            {
                case KTzParameters.netType:
                    val = Enum.Parse(typeof(KTzV2.Maths.Matrices.AdjacencyMatrix.AdjacencyMatrixType), value, true);
                    break;
                case KTzParameters.sType:
                    val = Enum.Parse(typeof(KTzV2.Synapses.SynapseType), value, true);
                    break;
                case KTzParameters.noiseType:
                    val = Enum.Parse(typeof(KTzV2.Synapses.NoiseType), value, true);
                    break;
                case KTzParameters.neuron:
                    val = Enum.Parse(typeof(KTzV2.Neurons.NeuronType), value, true);
                    break;
                case KTzParameters.stimType:
                    val = Enum.Parse(typeof(KTzV2.Stimuli.StimulusType), value, true);
                    break;
                case KTzParameters.iCond:
                    val = Enum.Parse(typeof(KTzV2.Data.InitialConditionType), value, true);
                    break;
                case KTzParameters.samp:
                    val = Enum.Parse(typeof(KTzV2.Data.NetworkSamplingType), value, true);
                    break;
                case KTzParameters.simType:
                    val = Enum.Parse(typeof(KTzV2.SimulationType), value, true);
                    break;
                case KTzParameters.dynType:
                    val = Enum.Parse(typeof(KTzV2.DynamicsSimType), value, true);
                    break;
                case KTzParameters.cVar:
                    val = Enum.Parse(typeof(KTzV2.Data.CountVariable), value, true);
                    break;
                case KTzParameters.bifWrite:
                    val = Enum.Parse(typeof(KTzV2.BifurcationWritePolicy), value, true);
                    break;
                case KTzParameters.indChoice:
                    val = Enum.Parse(typeof(KTzV2.Sims.Network.StimulusIndexChoice), value, true);
                    break;
                case KTzParameters.outAvgMode:
                    val = Enum.Parse(typeof(KTzV2.OutputAverageMode), value, true);
                    break;
                case KTzParameters.coupParam:
                    val = Enum.Parse(typeof(KTzV2.Sims.Network.CouplingParam), value, true);
                    break;
                case KTzParameters.oFileFormat:
                    val = Enum.Parse(typeof(KTzV2.Data.Header.OutputFileFormat), value, true);
                    break;
                case KTzParameters.ParamForRange:
                    val = Enum.Parse(typeof(KTzV2.ParamForRangeInDynamicsSim), value, true);
                    break;
                case KTzParameters.simTimeScheme:
                    val = Enum.Parse(typeof(KTzV2.Sims.Network.SimulationTimeScheme), value, true);
                    break;
                case KTzParameters.avgInp:
                case KTzParameters.wCSV:
                case KTzParameters.wData:
                case KTzParameters.wDif:
                case KTzParameters.wAvalDist:
                case KTzParameters.saveSpikeTimes:
                case KTzParameters.netDir:
                case KTzParameters.wObs:
                case KTzParameters.writeRhoTS:
                    val = Enum.Parse(typeof(KTzV2.Data.YesOrNoAnswer), value, true);
                    break;
                default:
                    if (KTzHeader.ParamList_String.ContainsKey(p))
                    {
                        if ((p == KTzParameters.netFile) || (p == KTzParameters.oFile))
                        {
                            char[] cseq;
                            if (value.Contains(System.IO.Path.DirectorySeparatorChar))
                                cseq = System.IO.Path.GetInvalidPathChars();
                            else
                                cseq = System.IO.Path.GetInvalidFileNameChars();
                            foreach (char c in cseq)
                            {
                                if (value.Contains(c))
                                {
                                    throw new ArgumentOutOfRangeException(String.Format("The value {0} has invalid characters.", value));
                                }
                            }
                        }
                        KTzHeader.SetPar(p, value);
                        return;
                    }
                    val = null;
                    break;
            }
            if (KTzHeader.IsDisorderInputParam(p))
            {
                KTzHeader.AddParDisorder(p, value);
            }
            else
            {
                if (val == null)
                {
                    throw new ArgumentOutOfRangeException(String.Format("Parameter {0} was not found as a text settable parameter.", p.ToString()));
                }
                if (KTzHeader.ParamList_Int32.ContainsKey(p))
                {
                    KTzHeader.SetPar(p, (Int32)val);
                    return;
                }
                throw new ArgumentOutOfRangeException(String.Format("Parameter {0} was not part of the Int32 valued list nor the String valued list.", p.ToString()));
            }
        }

        public static void AskForParamsChanges()
        {
            String input;
            KTzParameters p;
            Object val = new Object();

            Console.WriteLine("-");
            Console.WriteLine("to change any parameter, type its name and press <ENTER>");
            Console.WriteLine("or just press <ENTER> to continue");
            do
            {
                Console.Write("... param name: ");
                input = Console.ReadLine();
                try
                {
                    p = (KTzParameters)Enum.Parse(typeof(KTzParameters), input, false);
                    if (ParamList_Double.ContainsKey(p))
                    {
                        val = KTzHeader.query<Double>(input, ParamList_Double[p]);
                        KTzHeader.SetPar(p, (Double)val);
                    }
                    else
                    {
                        val = KTzHeader.query<Int32>(input, ParamList_Int32[p]);
                        KTzHeader.SetPar(p, (Int32)val);
                    }
                }
                catch (ArgumentException)
                {
                    KTzHeader.ErrorUnknownParameter(input, "");
                }
                Console.WriteLine("param set -> " + input + " = " + (String)val);
            } while (input != String.Empty);
        }

        public static T query<T>(String queryTxt, T defaultInput) where T : struct
        {
            System.ComponentModel.TypeConverter converter = System.ComponentModel.TypeDescriptor.GetConverter(typeof(T));
            if (!converter.CanConvertFrom(typeof(String)))
                throw new NotSupportedException("Can't convert to " + typeof(T).Name);

            String input;
            T result;

            while (true)
            {
                Console.WriteLine(" {0} (press only <ENTER> for default)", queryTxt);
                Console.Write(" default value = {0}, new value = ", defaultInput);
                input = Console.ReadLine();
                if (input != String.Empty)
                {
                    try
                    {
                        result = (T)converter.ConvertFromInvariantString(input);
                    }
                    catch (Exception e)
                    {
                        Console.WriteLine(e.Message);
                        continue;
                    }
                }
                else
                {
                    result = defaultInput;
                    Console.WriteLine("using default value: " + defaultInput);
                }
                break;
            }
            return result;
        }

        private static void AddParDisorder(KTzParameters parDisorder, String disorder_setting)
        {
            var pp = (KTzParameters)Enum.Parse(typeof(KTzParameters), parDisorder.ToString().Replace("Disorder", ""), false);
            KTzHeader.ParamDisorderSetting.Add(pp, new QuenchedDisorderParse(pp, disorder_setting));
            KTzHeader.AddPar(parDisorder, KTzHeader.ParamDisorderSetting[pp].ToString(), KTzHeader.ParamGroups[pp]);
            KTzHeader.ParamDescription.Add(parDisorder, String.Format("Quenched disorder definition for {0} parameter {1}", KTzHeader.ParamGroups[pp], pp));
        }

        private static void AddPar(KTzParameters par, Double val, KTzParamGroups group, bool allowsDisorder = false)
        {
            KTzHeader.ParamList_Double.Add(par, val);
            KTzHeader.ParamGroups.Add(par, group);
            if (allowsDisorder)
                KTzHeader.ParamAllowedQuenchedDisorder.Add(par);
        }

        private static void AddPar(KTzParameters par, Int32 val, KTzParamGroups group)
        {
            KTzHeader.ParamList_Int32.Add(par, val);
            KTzHeader.ParamGroups.Add(par, group);
        }

        private static void AddPar(KTzParameters par, Int32 val, KTzParamGroups group, Type T)
        {
            KTzHeader.ParamList_Int32.Add(par, val);
            KTzHeader.ParamGroups.Add(par, group);
            if (T.IsEnum)
                KTzHeader.ParamConverter.Add(par,k => Enum.GetName(T, k));
        }

        private static void AddPar(KTzParameters par, String val, KTzParamGroups group)
        {
            KTzHeader.ParamList_String.Add(par, val);
            KTzHeader.ParamGroups.Add(par, group);
        }

        public static void SetPar(KTzParameters par, Double val)
        {
            KTzHeader.ParamList_Double[par] = val;
            if (KTzHeader.ParamDisorderSetting.Keys.Contains(par))
                KTzHeader.ParamDisorderSetting[par].Set(val);
        }

        public static void SetPar(KTzParameters par, Int32 val)
        {
            KTzHeader.ParamList_Int32[par] = val;
        }

        public static void SetPar(KTzParameters par, String val)
        {
            KTzHeader.ParamList_String[par] = val;
        }
        
        /*
        public static T GetPar<T>(KTzParameters p)
        {
            System.ComponentModel.TypeConverter converter = System.ComponentModel.TypeDescriptor.GetConverter(typeof(T));
            T res;
            switch (Type.GetTypeCode(typeof(T)))
            {
                case TypeCode.Int32:
                    res = (T)converter.ConvertFrom(GetPar_Int32(p));
                    break;
                case TypeCode.Double:
                    res = (T)converter.ConvertFrom(GetPar_Double(p));
                    break;
                case TypeCode.String:
                    res = (T)converter.ConvertFrom(GetPar_String(p));
                    break;
                default:
                    throw new ArgumentOutOfRangeException("T","unknown type T");
            }
            return res;
        }
        /**/

        public static Double GetPar_Double(KTzParameters p)
        {
            if (KTzHeader.ParamDisorderSetting.Keys.Contains(p))
                return KTzHeader.ParamDisorderSetting[p].NextSample();
            return KTzHeader.ParamList_Double[p];
        }

        public static Int32 GetPar_Int32(KTzParameters p)
        {
            return KTzHeader.ParamList_Int32[p];
        }

        public static String GetPar_String(KTzParameters p)
        {
            return KTzHeader.ParamList_String[p];
        }

        public static Boolean HasParam(KTzParameters p)
        {
            return KTzHeader.ParamList_Double.ContainsKey(p) || KTzHeader.ParamList_Int32.ContainsKey(p) || KTzHeader.ParamList_String.ContainsKey(p);
        }

        private static void ErrorUnacceptedValue(String par, String val)
        {
            Console.WriteLine("Parameter {0} does not accept value {1}", par, val);
            Console.WriteLine("Exiting...");
            KTzHeader.ShowHelp();
        }

        private static void ErrorUnknownParameter(String par, String val)
        {
            Console.WriteLine("Unknown parameter's value: {0} = {1}", par, val);
            Console.WriteLine("Exiting...");
            KTzHeader.ShowHelp();
        }

        private static void ShowHelp()
        {
            Console.WriteLine("This program may run two different kinds of simulation: (1) Avalanche Dynamics and (2) Bifurcation.");
            Console.WriteLine("(1) The avalanche dynamics runs the specified network and returns avalanche distributions.");
            Console.WriteLine("(2) Bifurcation simulation runs the network for the specified J (or alpha) and I (or r) intervals and returns phase transition information as the parameters change.");
            Console.WriteLine("-");
            Console.WriteLine("USAGE");
            Console.WriteLine("-----");
            Console.WriteLine("$ " + KTzHeader.GetKTzV2ExeName() + " [-run] [-ask] [-wait] [PARAM1=VALUE1 PARAM2=VALUE2 ...]");
            Console.WriteLine("arguments between braces are optional");
            Console.WriteLine("-run\t\truns the program; if it is not specified, this message is printed");
            Console.WriteLine("-ask\t\tasks for params change");
            Console.WriteLine("-wait\t\twait for ENTER on the end of the program");
            Console.WriteLine("-");
            Console.WriteLine("QUENCHED DISORDER");
            Console.WriteLine("-------- --------");
            Console.WriteLine("Some parameters allow quenched disorder (see in the description of each parameter below).");
            Console.WriteLine("To set disorder, use: XDisorder=Type:X_Min_OR_X_Mean:X_Max_OR_X_Stddev");
            Console.WriteLine("E.g.,");
            Console.WriteLine("$ " + KTzHeader.GetKTzV2ExeName() + " [...] XDisorder=Type:X_Min_OR_X_Mean:X_Max_OR_X_Stddev [...]");
            Console.WriteLine("where X is the exact name of the parameter that accepts disorder; Type is either Gaussian or Uniform;");
            Console.WriteLine("if Type is Uniform, then provide min and max values; if type is Gaussian, then provide mean and stddev values");
            Console.WriteLine("For example: you can put in the param list: dDisorder=Gaussian:0.001:0.0001 TDisorder=Uniform:0.2:0.5");
            Console.WriteLine("to add Gaussian quenched disorder to the neuron's d parameter with mean=0.001 and stddev=0.0001;");
            Console.WriteLine("and uniform quenched disorder to T between 0.2 and 0.5.");
            Console.WriteLine("NOTE: Gaussian disorder is NOT 'restricted' and may change the sign of the parameter for some neurons/synapses!");
            Console.WriteLine("      Use with caution specially in the synapse parameters!");
            Console.WriteLine("-");
            Console.WriteLine("ALLOWED PARAMETERS");
            Console.WriteLine("------- ----------");
            Console.WriteLine("*** parameters which are not specified will have the values below");
            Console.WriteLine("*** the following names are available for PARAMX:");
            Console.Write(KTzHeader.GetAllParamString(true));
            //Console.WriteLine("TEST");
        }

        public static Dictionary<KTzParamGroups, Dictionary<String, String[]>> GetAllParamPairsAsStr()
        {
            var pGroups = new Dictionary<KTzParamGroups, Dictionary<String, String[]>>();
            foreach (var group in Enum.GetValues(typeof(KTzParamGroups)))
            {
                pGroups.Add((KTzParamGroups)group, KTzHeader.GetAllParamPairsAsStr((KTzParamGroups)group));
            }
            return pGroups;
        }

        public static Dictionary<String, String[]> GetAllParamPairsAsStr(KTzParamGroups group)
        {
            Func<Int32, String> get_par;
            Dictionary<String, String[]> p_list = new Dictionary<String, String[]>();
            foreach (KeyValuePair<KTzParameters, Double> par in KTzHeader.ParamList_Double)
            {
                if (KTzHeader.ParamGroups[par.Key] == group)
                    p_list.Add(par.Key.ToString(), new String[] { par.Value.ToString(), KTzHeader.ParamDescription[par.Key] });
            }
            foreach (KeyValuePair<KTzParameters, Int32> par in KTzHeader.ParamList_Int32)
            {
                if (KTzHeader.ParamGroups[par.Key] == group)
                {
                    if (KTzHeader.ParamConverter.ContainsKey(par.Key))
                        get_par = KTzHeader.ParamConverter[par.Key];
                    else
                        get_par = k => k.ToString();
                    p_list.Add(par.Key.ToString(), new String[] { get_par(par.Value), KTzHeader.ParamDescription[par.Key] });
                }
            }
            foreach (KeyValuePair<KTzParameters, String> par in KTzHeader.ParamList_String)
            {
                if (KTzHeader.ParamGroups[par.Key] == group)
                    p_list.Add(par.Key.ToString(), new String[] { par.Value, KTzHeader.ParamDescription[par.Key] });
            }
            return p_list;
        }

        public static String GetAllParamString(bool addLinePadding = false)
        {
            String str = "# " + KTzHeader.inputArgsStr + Environment.NewLine + "# -" + Environment.NewLine;
            foreach (var pGroup in KTzHeader.GetAllParamPairsAsStr())
            {
                String dashes = new String('-', (pGroup.Key.ToString() + " Parameters  ").Length);
                str += "#" + Environment.NewLine;
                str += "#-----------------" + dashes + "------------------" + Environment.NewLine;
                str += "#-----------------" + dashes + "------------------" + Environment.NewLine;
                str += "#----------------- " + pGroup.Key.ToString() + " Parameters ------------------" + Environment.NewLine;
                str += "#-----------------" + dashes + "------------------" + Environment.NewLine;
                str += "#-----------------" + dashes + "------------------" + Environment.NewLine;
                foreach (KeyValuePair<String, String[]> par in pGroup.Value)
                {
                    if (addLinePadding)
                        str += "# " + Environment.NewLine;
                    str += "# " + par.Key + " = " + par.Value[0] + "\t\t" + par.Value[1] + Environment.NewLine;
                    if (addLinePadding)
                        str += "# " + Environment.NewLine;
                }
            }
            return str;
        }

        public static Dictionary<KTzParamGroups, MatFileHandler.IVariable> GetAllParamPairsAsMatlabStruct(MatFileHandler.DataBuilder matDataBuilder)
        {
            //var matDataBuilder = new MatFileHandler.DataBuilder();
            var pGroups = new Dictionary<KTzParamGroups, MatFileHandler.IVariable>();
            foreach (var group in Enum.GetValues(typeof(KTzParamGroups)))
            {
                var par_in_group = KTzHeader.GetAllParamPairsAsMatlabStruct((KTzParamGroups)group, matDataBuilder);
                var s = matDataBuilder.NewStructureArray(par_in_group.Keys.Select(k => KTzHeader.KeepLetterOrDigit(k.ToString())), 1, 1);
                foreach (KTzParameters f in par_in_group.Keys)
                {
                    s[f.ToString(), 0, 0] = par_in_group[f];
                }
                pGroups.Add((KTzParamGroups)group, matDataBuilder.NewVariable((String)group.ToString() + "_Param",s));
            }
            return pGroups;
        }

        public static Dictionary<KTzParameters, MatFileHandler.IArray> GetAllParamPairsAsMatlabStruct(KTzParamGroups group, MatFileHandler.DataBuilder matDataBuilder)
        {
            var p_list = new Dictionary<KTzParameters, MatFileHandler.IArray>();
            foreach (KeyValuePair<KTzParameters, Double> par in KTzHeader.ParamList_Double)
            {
                if (KTzHeader.ParamGroups[par.Key] == group)
                    p_list.Add(par.Key, matDataBuilder.NewArray(new Double[] { par.Value }, 1, 1));
            }
            foreach (KeyValuePair<KTzParameters, Int32> par in KTzHeader.ParamList_Int32)
            {
                if (KTzHeader.ParamGroups[par.Key] == group)
                {
                    if (KTzHeader.ParamConverter.ContainsKey(par.Key))
                        p_list.Add(par.Key, matDataBuilder.NewCharArray(KTzHeader.ParamConverter[par.Key](par.Value)));
                    else
                        p_list.Add(par.Key, matDataBuilder.NewArray(new Int32[] { par.Value }, 1, 1));
                }
            }
            foreach (KeyValuePair<KTzParameters, String> par in KTzHeader.ParamList_String)
            {
                if (KTzHeader.ParamGroups[par.Key] == group)
                    p_list.Add(par.Key, matDataBuilder.NewCharArray(par.Value));
            }
            return p_list;
        }

        public static void ResetAllParamLists()
        {
            KTzHeader.ParamList_Double = new Dictionary<KTzParameters, Double>();
            KTzHeader.ParamList_Int32 = new Dictionary<KTzParameters, Int32>();
            KTzHeader.ParamList_String = new Dictionary<KTzParameters, String>();
            KTzHeader.ParamDescription = new Dictionary<KTzParameters, String>();
            KTzHeader.ParamGroups = new Dictionary<KTzParameters, KTzParamGroups>();
            KTzHeader.ParamConverter = new Dictionary<KTzParameters, Func<Int32, String>>();
            KTzHeader.ParamAllowedQuenchedDisorder = new HashSet<KTzParameters>();
            KTzHeader.ParamDisorderSetting = new Dictionary<KTzParameters, QuenchedDisorderParse>();

            KTzHeader.ParamDescription.Add(KTzParameters.K, "(allows quenched disorder) K self-interaction intensity between y and x");
            KTzHeader.AddPar(KTzParameters.K, 0.6, KTzParamGroups.Neuron, allowsDisorder: true);

            KTzHeader.ParamDescription.Add(KTzParameters.T, "(allows quenched disorder) T neuronal gain of the x variable");
            KTzHeader.AddPar(KTzParameters.T, 0.35, KTzParamGroups.Neuron, allowsDisorder: true);

            KTzHeader.ParamDescription.Add(KTzParameters.d, "(allows quenched disorder) delta: z recovery inverse time scale, controls refractory period and burst damping");
            KTzHeader.AddPar(KTzParameters.d, 0.001, KTzParamGroups.Neuron, allowsDisorder: true);

            KTzHeader.ParamDescription.Add(KTzParameters.l, "(allows quenched disorder) lambda: z-x coupling time scale, controls refractory period");
            KTzHeader.AddPar(KTzParameters.l, 0.008, KTzParamGroups.Neuron, allowsDisorder: true);

            KTzHeader.ParamDescription.Add(KTzParameters.xR, "(allows quenched disorder) z-x coupling recovery x potential, controls burst duration");
            KTzHeader.AddPar(KTzParameters.xR, -0.7, KTzParamGroups.Neuron, allowsDisorder: true);

            KTzHeader.ParamDescription.Add(KTzParameters.H, "(allows quenched disorder) polarizing current (external field in x; or in y if KTz2H neuron is selected)");
            KTzHeader.AddPar(KTzParameters.H, 0.0D, KTzParamGroups.Neuron, allowsDisorder: true);

            KTzHeader.ParamDescription.Add(KTzParameters.Q, "(allows quenched disorder) polarizing current for x in KTz2Tanh (external field)");
            KTzHeader.AddPar(KTzParameters.Q, 0.0D, KTzParamGroups.Neuron, allowsDisorder: true);

            KTzHeader.ParamDescription.Add(KTzParameters.J, "(allows quenched disorder) coupling intensity (conductance)");
            KTzHeader.AddPar(KTzParameters.J, -0.15, KTzParamGroups.Synapse, allowsDisorder: true);

            KTzHeader.ParamDescription.Add(KTzParameters.tauf, "(allows quenched disorder) tau_f: recovery time scale of the f variable for the chemical synapse");
            KTzHeader.AddPar(KTzParameters.tauf, 2.0, KTzParamGroups.Synapse, allowsDisorder: true);

            KTzHeader.ParamDescription.Add(KTzParameters.taug, "(allows quenched disorder) tau_g: recovery time scale of the g variable for the chemical synapse");
            KTzHeader.AddPar(KTzParameters.taug, 2.0, KTzParamGroups.Synapse, allowsDisorder: true);

            KTzHeader.ParamDescription.Add(KTzParameters.x0, "neuron initial condition on x");
            KTzHeader.AddPar(KTzParameters.x0, -0.6971564118917724, KTzParamGroups.InitCond);

            KTzHeader.ParamDescription.Add(KTzParameters.y0, "neuron initial condition on y");
            KTzHeader.AddPar(KTzParameters.y0, -0.6971564118917724, KTzParamGroups.InitCond);

            KTzHeader.ParamDescription.Add(KTzParameters.z0, "neuron initial condition on z");
            KTzHeader.AddPar(KTzParameters.z0, -0.0227487048658225, KTzParamGroups.InitCond);

            KTzHeader.ParamDescription.Add(KTzParameters.minJ, "minimum J on bifurcation simulations");
            KTzHeader.AddPar(KTzParameters.minJ, -0.2, KTzParamGroups.Synapse);

            KTzHeader.ParamDescription.Add(KTzParameters.maxJ, "maximum J on bifurcation simulations");
            KTzHeader.AddPar(KTzParameters.maxJ, -0.01, KTzParamGroups.Synapse);

            KTzHeader.ParamDescription.Add(KTzParameters.I, "stimulus intensity");
            KTzHeader.AddPar(KTzParameters.I, 0.1, KTzParamGroups.Stimulus);

            KTzHeader.ParamDescription.Add(KTzParameters.minI, "minimum I on bifurcation simulations");
            KTzHeader.AddPar(KTzParameters.minI, 0.03, KTzParamGroups.Stimulus);

            KTzHeader.ParamDescription.Add(KTzParameters.maxI, "maximum I on bifurcation simulations");
            KTzHeader.AddPar(KTzParameters.maxI, 0.13, KTzParamGroups.Stimulus);

            KTzHeader.ParamDescription.Add(KTzParameters.r, "Poisson stimulus rate");
            KTzHeader.AddPar(KTzParameters.r, 0.1, KTzParamGroups.Stimulus);

            KTzHeader.ParamDescription.Add(KTzParameters.minr, "minimum Poisson stimulus rate on bifurcation simulations");
            KTzHeader.AddPar(KTzParameters.minr, 0.00001, KTzParamGroups.Stimulus);

            KTzHeader.ParamDescription.Add(KTzParameters.maxr, "maximum Poisson stimulus rate on bifurcation simulations");
            KTzHeader.AddPar(KTzParameters.maxr, 10.0, KTzParamGroups.Stimulus);

            KTzHeader.ParamDescription.Add(KTzParameters.IStdDev, "std deviation of I in case of Poisson stimulus: set to 0 to have the same intensity for every stimulus");
            KTzHeader.AddPar(KTzParameters.IStdDev, 0.0, KTzParamGroups.Stimulus);

            KTzHeader.ParamDescription.Add(KTzParameters.R, "(allows quenched disorder) noise amplitude");
            KTzHeader.AddPar(KTzParameters.R, 0.036, KTzParamGroups.Synapse, allowsDisorder: true);

            KTzHeader.ParamDescription.Add(KTzParameters.noiseRatio, "(allows quenched disorder) fraction of J to use as noise amplitude");
            KTzHeader.AddPar(KTzParameters.noiseRatio, 0.1, KTzParamGroups.Synapse, allowsDisorder: true);

            KTzHeader.ParamDescription.Add(KTzParameters.alpha, "(allows quenched disorder) LHG-like dynamics for J in Chemical Synapses: coupling intensity * u");
            KTzHeader.AddPar(KTzParameters.alpha, 1.5, KTzParamGroups.Synapse, allowsDisorder: true);

            KTzHeader.ParamDescription.Add(KTzParameters.minalpha, "LHG-like dynamics for J in Chemical Synapses: minimum couplingintensity*u on bifurcation simulations");
            KTzHeader.AddPar(KTzParameters.minalpha, 0.8, KTzParamGroups.Synapse);

            KTzHeader.ParamDescription.Add(KTzParameters.maxalpha, "LHG-like dynamics for J in Chemical Synapses: maximum couplingintensity*u on bifurcation simulations");
            KTzHeader.AddPar(KTzParameters.maxalpha, 1.8, KTzParamGroups.Synapse);

            KTzHeader.ParamDescription.Add(KTzParameters.u, "(allows quenched disorder) LHG-like dynamics for J in Chemical Synapses: synapse depression strength");
            KTzHeader.AddPar(KTzParameters.u, 0.2, KTzParamGroups.Synapse, allowsDisorder: true);

            KTzHeader.ParamDescription.Add(KTzParameters.tauJ, "(allows quenched disorder) LHG-like dynamics for J in Chemical Synapses: coupling J recovery time scale");
            KTzHeader.AddPar(KTzParameters.tauJ, 3000.0, KTzParamGroups.Synapse, allowsDisorder: true);

            KTzHeader.ParamDescription.Add(KTzParameters.dt, "LHG-like dynamics for J in Chemical Synapses: precision on the integration of the J equations (with Euler method)");
            KTzHeader.AddPar(KTzParameters.dt, 0.5, KTzParamGroups.Synapse);

            KTzHeader.ParamDescription.Add(KTzParameters.rewP, "Watts-Strogatz or Random Graph Network parameter - rewire probability or edge creation prob");
            KTzHeader.AddPar(KTzParameters.rewP, 0.02, KTzParamGroups.Network);

            KTzHeader.ParamDescription.Add(KTzParameters.sampFrac, "fraction of neurons to sample from the network -- when subsampling is intended");
            KTzHeader.AddPar(KTzParameters.sampFrac, 1.0, KTzParamGroups.Simulation);

            KTzHeader.ParamDescription.Add(KTzParameters.Lx, "Amount of elements on x axis (linear size) (total amount is always assumed as Lx*Ly*Lz)");
            KTzHeader.AddPar(KTzParameters.Lx, 40, KTzParamGroups.Network);

            KTzHeader.ParamDescription.Add(KTzParameters.Ly, "Amount of elements on y axis (linear size) (total amount is always assumed as Lx*Ly*Lz)");
            KTzHeader.AddPar(KTzParameters.Ly, 40, KTzParamGroups.Network);

            KTzHeader.ParamDescription.Add(KTzParameters.Lz, "Amount of elements on z axis (linear size) (total amount is always assumed as Lx*Ly*Lz)");
            KTzHeader.AddPar(KTzParameters.Lz, 1, KTzParamGroups.Network);

            KTzHeader.ParamDescription.Add(KTzParameters.nSteps, "total timesteps to run each simulation");
            KTzHeader.AddPar(KTzParameters.nSteps, 10000, KTzParamGroups.Simulation);

            KTzHeader.ParamDescription.Add(KTzParameters.nStart, "timestep in which the network starts to be recorded");
            KTzHeader.AddPar(KTzParameters.nStart, 0, KTzParamGroups.Simulation);

            KTzHeader.ParamDescription.Add(KTzParameters.netType, "network architecture" + Environment.NewLine + KTzHeader.GetEnumListStr<KTzV2.Maths.Matrices.AdjacencyMatrix.AdjacencyMatrixType>());
            KTzHeader.AddPar(KTzParameters.netType, (Int32)KTzV2.Maths.Matrices.AdjacencyMatrix.AdjacencyMatrixType.SquareLatticeFreeBC, KTzParamGroups.Network, typeof(KTzV2.Maths.Matrices.AdjacencyMatrix.AdjacencyMatrixType));

            KTzHeader.ParamDescription.Add(KTzParameters.netDir, "directed = yes: the synapses are one-way i->j; directed = no: two-way synapses" + Environment.NewLine + KTzHeader.GetEnumListStr<KTzV2.Data.YesOrNoAnswer>());
            KTzHeader.AddPar(KTzParameters.netDir, (Int32)KTzV2.Data.YesOrNoAnswer.No, KTzParamGroups.Network, typeof(KTzV2.Data.YesOrNoAnswer));

            KTzHeader.ParamDescription.Add(KTzParameters.sType, "type of the synapse" + Environment.NewLine + KTzHeader.GetEnumListStr<KTzV2.Synapses.SynapseType>());
            KTzHeader.AddPar(KTzParameters.sType, (Int32)KTzV2.Synapses.SynapseType.KTNoisyChemicalSynapse, KTzParamGroups.Synapse, typeof(KTzV2.Synapses.SynapseType));

            KTzHeader.ParamDescription.Add(KTzParameters.noiseType, "type of the synaptic noise" + Environment.NewLine + KTzHeader.GetEnumListStr<KTzV2.Synapses.NoiseType>());
            KTzHeader.AddPar(KTzParameters.noiseType, (Int32)KTzV2.Synapses.NoiseType.ProportionalAmplitude, KTzParamGroups.Synapse, typeof(KTzV2.Synapses.NoiseType));

            KTzHeader.ParamDescription.Add(KTzParameters.neuron, "type of neuron to use" + Environment.NewLine + KTzHeader.GetEnumListStr<KTzV2.Neurons.NeuronType>());
            KTzHeader.AddPar(KTzParameters.neuron, (Int32)KTzV2.Neurons.NeuronType.KTz, KTzParamGroups.Neuron, typeof(KTzV2.Neurons.NeuronType));

            KTzHeader.ParamDescription.Add(KTzParameters.iCond, "type of initial condition" + Environment.NewLine + KTzHeader.GetEnumListStr<KTzV2.Data.InitialConditionType>());
            KTzHeader.AddPar(KTzParameters.iCond, (Int32)KTzV2.Data.InitialConditionType.ProgramSpecified, KTzParamGroups.InitCond, typeof(KTzV2.Data.InitialConditionType));

            KTzHeader.ParamDescription.Add(KTzParameters.dim, "dimension of the network (network will have L^dim neurons)");
            KTzHeader.AddPar(KTzParameters.dim, 2, KTzParamGroups.Network);

            KTzHeader.ParamDescription.Add(KTzParameters.stimType, "type of the stimulus" + Environment.NewLine + KTzHeader.GetEnumListStr<KTzV2.Stimuli.StimulusType>());
            KTzHeader.AddPar(KTzParameters.stimType, (Int32)KTzV2.Stimuli.StimulusType.Delta, KTzParamGroups.Stimulus, typeof(KTzV2.Stimuli.StimulusType));

            KTzHeader.ParamDescription.Add(KTzParameters.simTimeScheme, "CAUTION: simulation will become VERY SLOW for small Poisson Rate r; if 'ProportionalToPoissonRate' then (nSteps - nStartStep) = 10/(r*N), where r-> Poisson rate; N=Lx*Ly*Lz (number of neurons)" + Environment.NewLine + KTzHeader.GetEnumListStr<KTzV2.Sims.Network.SimulationTimeScheme>());
            KTzHeader.AddPar(KTzParameters.simTimeScheme, (Int32)KTzV2.Sims.Network.SimulationTimeScheme.Free, KTzParamGroups.Simulation, typeof(KTzV2.Sims.Network.SimulationTimeScheme));

            KTzHeader.ParamDescription.Add(KTzParameters.sStim, "timestep in which the network will be stimulated for each run");
            KTzHeader.AddPar(KTzParameters.sStim, 0, KTzParamGroups.Stimulus);

            KTzHeader.ParamDescription.Add(KTzParameters.iStim, "index of the neuron which will be stimulated");
            KTzHeader.AddPar(KTzParameters.iStim, -1, KTzParamGroups.Stimulus);

            KTzHeader.ParamDescription.Add(KTzParameters.tBin, "size of the bin (in timesteps) to divide nSteps in order to run the simulation");
            KTzHeader.AddPar(KTzParameters.tBin, 20, KTzParamGroups.Simulation);

            KTzHeader.ParamDescription.Add(KTzParameters.wData, "choose whether the program will run a simulation only to write output data in the case of dynamics simulation" + Environment.NewLine + KTzHeader.GetEnumListStr<KTzV2.Data.YesOrNoAnswer>());
            KTzHeader.AddPar(KTzParameters.wData, (Int32)KTzV2.Data.YesOrNoAnswer.No, KTzParamGroups.Simulation, typeof(KTzV2.Data.YesOrNoAnswer));

            KTzHeader.ParamDescription.Add(KTzParameters.wCSV, "if wData=Yes, then wCSV=Yes writes a CSV file with the x_i(t) data" + Environment.NewLine + KTzHeader.GetEnumListStr<KTzV2.Data.YesOrNoAnswer>());
            KTzHeader.AddPar(KTzParameters.wCSV, (Int32)KTzV2.Data.YesOrNoAnswer.No, KTzParamGroups.Simulation, typeof(KTzV2.Data.YesOrNoAnswer));

            KTzHeader.ParamDescription.Add(KTzParameters.wAvalDist, "choose whether the program will run a simulation only to write an avalanche spike distribution file in the case of dynamics simulation" + Environment.NewLine + KTzHeader.GetEnumListStr<KTzV2.Data.YesOrNoAnswer>());
            KTzHeader.AddPar(KTzParameters.wAvalDist, (Int32)KTzV2.Data.YesOrNoAnswer.Yes, KTzParamGroups.Simulation, typeof(KTzV2.Data.YesOrNoAnswer));

            KTzHeader.ParamDescription.Add(KTzParameters.saveSpikeTimes, "chooses whether to save spike times for each neuron during a simulation" + Environment.NewLine + KTzHeader.GetEnumListStr<KTzV2.Data.YesOrNoAnswer>());
            KTzHeader.AddPar(KTzParameters.saveSpikeTimes, (Int32)KTzV2.Data.YesOrNoAnswer.Yes, KTzParamGroups.Simulation, typeof(KTzV2.Data.YesOrNoAnswer));

            KTzHeader.ParamDescription.Add(KTzParameters.wObs, "write observations of the avalanche sizes for Bifurcation sim type for each (par1,par2) pair" + Environment.NewLine + KTzHeader.GetEnumListStr<KTzV2.Data.YesOrNoAnswer>());
            KTzHeader.AddPar(KTzParameters.wObs, (Int32)KTzV2.Data.YesOrNoAnswer.No, KTzParamGroups.Simulation, typeof(KTzV2.Data.YesOrNoAnswer));

            KTzHeader.ParamDescription.Add(KTzParameters.nJ, "amount of J (synaptic coupling) on bifurcation simulations");
            KTzHeader.AddPar(KTzParameters.nJ, 100, KTzParamGroups.Synapse);

            KTzHeader.ParamDescription.Add(KTzParameters.nI, "amount of I (stimulus intensity) on bifurcation simulations");
            KTzHeader.AddPar(KTzParameters.nI, 100, KTzParamGroups.Stimulus);

            KTzHeader.ParamDescription.Add(KTzParameters.nr, "amount of r (Poisson rate) on bifurcation simulations");
            KTzHeader.AddPar(KTzParameters.nr, 100, KTzParamGroups.Stimulus);

            KTzHeader.ParamDescription.Add(KTzParameters.nalpha, "amount of alpha (LHG J adaptation intensity) on bifurcation simulations");
            KTzHeader.AddPar(KTzParameters.nalpha, 100, KTzParamGroups.Synapse);

            KTzHeader.ParamDescription.Add(KTzParameters.wDif, "choose whether the program will write a file containing xi-xj data; WARNING FILE WILL BE VERY LARGE if N is big" + Environment.NewLine + KTzHeader.GetEnumListStr<KTzV2.Data.YesOrNoAnswer>());
            KTzHeader.AddPar(KTzParameters.wDif, (Int32)KTzV2.Data.YesOrNoAnswer.No, KTzParamGroups.Simulation, typeof(KTzV2.Data.YesOrNoAnswer));

            KTzHeader.ParamDescription.Add(KTzParameters.nSim, "amount of realizations for each pair (stimulus,coupling) on bifurcation simulations");
            KTzHeader.AddPar(KTzParameters.nSim, 20, KTzParamGroups.Simulation);

            KTzHeader.ParamDescription.Add(KTzParameters.deltaT, "time interval between two consecutive delta stimuli");
            KTzHeader.AddPar(KTzParameters.deltaT, 20, KTzParamGroups.Stimulus);

            KTzHeader.ParamDescription.Add(KTzParameters.nNeigh, "Watts-Strogatz network parameter - amount of neighbours in the initial configuration of WS network");
            KTzHeader.AddPar(KTzParameters.nNeigh, 4, KTzParamGroups.Network);

            KTzHeader.ParamDescription.Add(KTzParameters.nConn, "Barabasi-Albert network parameter - amount of connections which a new node will initially have when creating BA network");
            KTzHeader.AddPar(KTzParameters.nConn, 3, KTzParamGroups.Network);

            KTzHeader.ParamDescription.Add(KTzParameters.cVar, "choose if the program will count number of neurons spiking or number of spikes; if NumberOfNeurons is set, then each of the nSim run is performed during nSteps timesteps" + Environment.NewLine + KTzHeader.GetEnumListStr<KTzV2.Data.CountVariable>());
            KTzHeader.AddPar(KTzParameters.cVar, (Int32)KTzV2.Data.CountVariable.NumberOfNeurons, KTzParamGroups.Simulation, typeof(KTzV2.Data.CountVariable));

            KTzHeader.ParamDescription.Add(KTzParameters.avgInp, "choose whether the neuron will average its synaptic input (may be useful in mean-field networks)" + Environment.NewLine + KTzHeader.GetEnumListStr<KTzV2.Data.YesOrNoAnswer>());
            KTzHeader.AddPar(KTzParameters.avgInp, (Int32)KTzV2.Data.YesOrNoAnswer.No, KTzParamGroups.Neuron, typeof(KTzV2.Data.YesOrNoAnswer));

            KTzHeader.ParamDescription.Add(KTzParameters.samp, "choose the sampling type (full or subsampled)" + Environment.NewLine + KTzHeader.GetEnumListStr<KTzV2.Data.NetworkSamplingType>());
            KTzHeader.AddPar(KTzParameters.samp, (Int32)KTzV2.Data.NetworkSamplingType.Full, KTzParamGroups.Simulation, typeof(KTzV2.Data.NetworkSamplingType));

            KTzHeader.ParamDescription.Add(KTzParameters.rest, "interval (in tBin's) between an avalanche and a new stimulus for no activity stimulus type (i.e., if no activity is detected in a time bin)");
            KTzHeader.AddPar(KTzParameters.rest, 10, KTzParamGroups.Simulation);

            KTzHeader.ParamDescription.Add(KTzParameters.simType, "choose if the program will run bifurcation or dynamics" + Environment.NewLine + KTzHeader.GetEnumListStr<KTzV2.SimulationType>());
            KTzHeader.AddPar(KTzParameters.simType, (Int32)KTzV2.SimulationType.Dynamics, KTzParamGroups.Simulation, typeof(KTzV2.SimulationType));

            KTzHeader.ParamDescription.Add(KTzParameters.dynType, "each avalanche is generated independently of each other or after a time window" + Environment.NewLine + KTzHeader.GetEnumListStr<KTzV2.DynamicsSimType>());
            KTzHeader.AddPar(KTzParameters.dynType, (Int32)KTzV2.DynamicsSimType.NetworkReset, KTzParamGroups.Simulation, typeof(KTzV2.DynamicsSimType));

            KTzHeader.ParamDescription.Add(KTzParameters.oFile, "prefix of the output file name -- do not use an extension");
            KTzHeader.AddPar(KTzParameters.oFile, "", KTzParamGroups.Simulation);

            KTzHeader.ParamDescription.Add(KTzParameters.oFileFormat, "Chooses whether to write the main outputs either in txt (.dat extension) or mat-file format (.mat MATLAB format)" + Environment.NewLine + KTzHeader.GetEnumListStr<KTzV2.Data.Header.OutputFileFormat>());
            KTzHeader.AddPar(KTzParameters.oFileFormat, (Int32)KTzV2.Data.Header.OutputFileFormat.txt, KTzParamGroups.Simulation, typeof(KTzV2.Data.Header.OutputFileFormat));

            KTzHeader.ParamDescription.Add(KTzParameters.netFile, "input file with adjacency matrix");
            KTzHeader.AddPar(KTzParameters.netFile, "", KTzParamGroups.Simulation);

            KTzHeader.ParamDescription.Add(KTzParameters.JRange, "range values for J parameter when running Bifurcation simulation; if set ignores minJ, maxJ, nJ (comma separated list of: min:max:dx OR range(min:max:dx) OR linspace(min:max:nx) OR logspace(min:max:nx))");
            KTzHeader.AddPar(KTzParameters.JRange, "", KTzParamGroups.Synapse);

            KTzHeader.ParamDescription.Add(KTzParameters.IRange, "range values for I parameter when running Bifurcation simulation; if set ignores minI, maxI, nI (comma separated list of: min:max:dx OR range(min:max:dx) OR linspace(min:max:nx) OR logspace(min:max:nx))");
            KTzHeader.AddPar(KTzParameters.IRange, "", KTzParamGroups.Stimulus);

            KTzHeader.ParamDescription.Add(KTzParameters.rRange, "log-spaced range values for r parameter when running Bifurcation simulation; if set ignores minr, maxr, nr (comma separated list of: min:max:dx OR range(min:max:dx) OR linspace(min:max:nx) OR logspace(min:max:nx))");
            KTzHeader.AddPar(KTzParameters.rRange, "", KTzParamGroups.Stimulus);

            KTzHeader.ParamDescription.Add(KTzParameters.alphaRange, "range values for alpha parameter when running Bifurcation simulation; if set ignores minalpha, maxalpha, nalpha (comma separated list of: min:max:dx OR range(min:max:dx) OR linspace(min:max:nx) OR logspace(min:max:nx))");
            KTzHeader.AddPar(KTzParameters.alphaRange, "", KTzParamGroups.Synapse);
            
            KTzHeader.ParamDescription.Add(KTzParameters.ParamForRange, "parameter used to create a range for simulation of the type DynamicsWithinParamRange" + Environment.NewLine + KTzHeader.GetEnumListStr<KTzV2.ParamForRangeInDynamicsSim>());
            KTzHeader.AddPar(KTzParameters.ParamForRange, (Int32)KTzV2.ParamForRangeInDynamicsSim.J, KTzParamGroups.Simulation, typeof(KTzV2.ParamForRangeInDynamicsSim));

            KTzHeader.ParamDescription.Add(KTzParameters.bifWrite, "OnTheFly == writes data after each pair (I,J); OnTheEnd == writes after everything is finished" + Environment.NewLine + KTzHeader.GetEnumListStr<KTzV2.BifurcationWritePolicy>());
            KTzHeader.AddPar(KTzParameters.bifWrite, (Int32)KTzV2.BifurcationWritePolicy.InTheEnd, KTzParamGroups.Simulation, typeof(KTzV2.BifurcationWritePolicy));

            KTzHeader.ParamDescription.Add(KTzParameters.indChoice, "iStim!=-1 forces indChoice==Fixed; method for choosing stimulus index" + Environment.NewLine + KTzHeader.GetEnumListStr<KTzV2.Sims.Network.StimulusIndexChoice>());
            KTzHeader.AddPar(KTzParameters.indChoice, (Int32)KTzV2.Sims.Network.StimulusIndexChoice.Fixed, KTzParamGroups.Stimulus, typeof(KTzV2.Sims.Network.StimulusIndexChoice));

            KTzHeader.ParamDescription.Add(KTzParameters.outAvgMode, "determines whether the output average should be taken over many avalanches or just temporal average (for bifurcation simulations)" + Environment.NewLine + KTzHeader.GetEnumListStr<KTzV2.OutputAverageMode>());
            KTzHeader.AddPar(KTzParameters.outAvgMode, (Int32)KTzV2.OutputAverageMode.OverRealizations, KTzParamGroups.Simulation, typeof(KTzV2.OutputAverageMode));

            KTzHeader.ParamDescription.Add(KTzParameters.Theta, "threshold of SI(RS) element");
            KTzHeader.AddPar(KTzParameters.Theta, 0.1, KTzParamGroups.Neuron, allowsDisorder: true);

            KTzHeader.ParamDescription.Add(KTzParameters.coupParam, "determines whether the J are given by J parameter or by adjacency matrix" + Environment.NewLine + KTzHeader.GetEnumListStr<KTzV2.Sims.Network.CouplingParam>());
            KTzHeader.AddPar(KTzParameters.coupParam, (Int32)KTzV2.Sims.Network.CouplingParam.Homogeneous, KTzParamGroups.Synapse, typeof(KTzV2.Sims.Network.CouplingParam));

            KTzHeader.ParamDescription.Add(KTzParameters.ThetaRange, "range of values for threshold of SI(RS) element in bifurcation simulation  (comma separated list of: min:max:dx OR range(min:max:dx) OR linspace(min:max:nx) OR logspace(min:max:nx))");
            KTzHeader.AddPar(KTzParameters.ThetaRange, "", KTzParamGroups.Neuron);

            KTzHeader.ParamDescription.Add(KTzParameters.minTheta, "min value for threshold of SI(RS) element in bifurcation simulation");
            KTzHeader.AddPar(KTzParameters.minTheta, 0.0, KTzParamGroups.Neuron);

            KTzHeader.ParamDescription.Add(KTzParameters.maxTheta, "max value for threshold of SI(RS) element in bifurcation simulation");
            KTzHeader.AddPar(KTzParameters.maxTheta, 1.0, KTzParamGroups.Neuron);

            KTzHeader.ParamDescription.Add(KTzParameters.nTheta, "number of values for threshold of SI(RS) element in bifurcation simulation");
            KTzHeader.AddPar(KTzParameters.nTheta, 100, KTzParamGroups.Neuron);

            KTzHeader.ParamDescription.Add(KTzParameters.writeRhoTS, "writes rho(t) for each J in JRange (only if outAvgMode=OverTime); rho(t) = sum x(t)>0" + Environment.NewLine + KTzHeader.GetEnumListStr<KTzV2.Data.YesOrNoAnswer>());
            KTzHeader.AddPar(KTzParameters.writeRhoTS, (Int32)KTzV2.Data.YesOrNoAnswer.No, KTzParamGroups.Simulation, typeof(KTzV2.Data.YesOrNoAnswer));
        }

        public static void GetMinMaxNFromStr(String rangeStr, out Double xMin, out Double xMax, out Int32 nx)
        {
            bool contains_range = rangeStr.Contains("range");
            Double dx;
            String[] ssp = rangeStr.Replace("range(", "").Replace(")", "").Split(':');
            xMin = Convert.ToDouble(ssp[0]);
            if (ssp.Length == 2)
            {
                dx = 1.0D;
                xMax = Convert.ToDouble(ssp[1]);
            }
            else if (ssp.Length == 3)
            {
                if (contains_range)
                {
                    dx = Convert.ToDouble(ssp[2]);
                    xMax = Convert.ToDouble(ssp[1]);
                }
                else 
                {
                    dx = Convert.ToDouble(ssp[1]);
                    xMax = Convert.ToDouble(ssp[2]);
                }
            }
            else
            {
                dx = 1.0D;
                xMax = xMin;
            }
            nx = Convert.ToInt32(Math.Round(1.0D + (xMax - xMin) / dx, MidpointRounding.ToEven));
        }

        public static void GetMinMaxNFromLinspaceStr(String linspaceStr, out Double xMin, out Double xMax, out Int32 nx)
        {
            String[] ss_parts = linspaceStr.Replace("logspace(","").Replace("linspace(", "").Replace(")", "").Split(':');
            Func<String,String> get_number = s => System.Text.RegularExpressions.Regex.Replace(s, "[^.efEF0-9+-]", "");
            xMin = Convert.ToDouble(get_number(ss_parts[0]));
            xMax = Convert.ToDouble(get_number(ss_parts[1]));
            nx   = Convert.ToInt32( get_number(ss_parts[2]));
        }

        public static Double[] GetRangeFor(KTzParameters p)
        {
            if (!((p == KTzParameters.J) ||
                (p == KTzParameters.alpha) ||
                (p == KTzParameters.r) ||
                (p == KTzParameters.I) ||
                (p == KTzParameters.Theta)))
            {
                throw new ArgumentException(String.Format("KTzV2.Header: parameter {0} does not have a range!", p.ToString()));
            }
            String parName = p.ToString() + "Range";
            String s = KTzHeader.GetPar_String((KTzParameters)Enum.Parse(typeof(KTzParameters), parName));
            if (s == String.Empty)
            {
                return KTzHeader.GetMinMaxRange(p);
            }
            else
            {
                Double xMin, xMax;
                Int32 nx;
                List<Double> range = new List<Double>();
                String[] sr = s.Split(new char[] { ',' }, StringSplitOptions.RemoveEmptyEntries);
                foreach (String ss in sr)
                {
                    String myss = ss.ToLower();
                    if (myss.Contains("range"))
                    {
                        KTzHeader.GetMinMaxNFromStr(myss, out xMin, out xMax, out nx);
                        range.AddRange((p == KTzParameters.r ? KTzHeader.GetMinMaxLog10Range(xMin, xMax, nx) : KTzHeader.GetMinMaxRange(xMin, xMax, nx)));
                    }
                    else
                    {
                        if (myss.Contains("linspace"))
                        {
                            GetMinMaxNFromLinspaceStr(myss, out xMin, out xMax, out nx);
                            range.AddRange(KTzHeader.GetMinMaxRange(xMin, xMax, nx));
                        }
                        else if (myss.Contains("logspace"))
                        {
                            GetMinMaxNFromLinspaceStr(myss, out xMin, out xMax, out nx);
                            range.AddRange(KTzHeader.GetMinMaxLog10Range(xMin, xMax, nx));
                        }
                        else
                        {
                            if (myss.Contains(':')) // matlab-like range
                            {
                                KTzHeader.GetMinMaxNFromStr(myss, out xMin, out xMax, out nx);
                                range.AddRange((p == KTzParameters.r ? KTzHeader.GetMinMaxLog10Range(xMin, xMax, nx) : KTzHeader.GetMinMaxRange(xMin, xMax, nx)));
                            }
                            else
                                range.Add(Convert.ToDouble(myss));
                        }
                        
                    }
                }
                range.Sort();
                return range.ToArray();
            }
        }

        public static Double[] GetMinMaxLog10Range(Double pMin, Double pMax, Int32 np, Double nbase = 10.0D)
        {
            pMin = Math.Log(pMin) / Math.Log(nbase);
            pMax = Math.Log(pMax) / Math.Log(nbase);
            Double[] y_list = KTzHeader.GetMinMaxRange(pMin, pMax, np);
            int i = -1;
            foreach (Double y in y_list)
            {
                i++;
                y_list[i] = Math.Pow(nbase, y);
            }
            return y_list;
        }

        public static Double[] GetMinMaxRange(Double pMin, Double pMax, Int32 np)
        {
            Double[] range = new Double[np]; // +1 to include the Jmax as well
            Double dp = (np > 1 ? (pMax - pMin) / Convert.ToDouble(np - 1) : 0.0D);
            Int32 i = 0;
            while (i < np)
            {
                range[i] = pMin + (Double)i * dp;
                i++;
            }
            return range;
        }

        private static Double[] GetMinMaxRange(KTzParameters p)
        {
            Double pMin = KTzHeader.GetPar_Double((KTzParameters)Enum.Parse(typeof(KTzParameters), "min" + p.ToString()));
            Double pMax = KTzHeader.GetPar_Double((KTzParameters)Enum.Parse(typeof(KTzParameters), "max" + p.ToString()));
            Int32 np = KTzHeader.GetPar_Int32((KTzParameters)Enum.Parse(typeof(KTzParameters), "n" + p.ToString()));
            if (pMin == pMax)
                np = 1;
            if (p == KTzParameters.r)
            {
                return KTzHeader.GetMinMaxLog10Range(pMin, pMax, np);
            }
            else
            {
                return KTzHeader.GetMinMaxRange(pMin, pMax, np);
            }
        }

        public static String GetEnumListStr<T>()
        {
            return String.Join(
                       Environment.NewLine,
                       Enum.GetNames(typeof(T))
                           .ToList()
                           .Select(str =>
                               {
                                   var at = (T)Enum.Parse(typeof(T), str);
                                   return "#\t\t\t\t" + Convert.ToInt32(at).ToString() + " = " + str;
                               })
                            .ToArray<String>());
        }

        public static String CheckAndGetFileName(String fileNameWithExtension)
        {
            String[] temp = fileNameWithExtension.Split(new char[] { '.' });
            String fileName = String.Empty, fileExt;
            Int32 n = temp.Length - 1;
            Int32 i = 0;
            while (i < n)
            {
                fileName += temp[i] + ".";
                i++;
            }
            fileExt = temp[n];
            return KTzHeader.CheckAndGetFileName(fileName.Substring(0, fileName.Length - 1), fileExt);
        }

        public static String CheckAndGetFileName(String fileName, String ext)
        {
            // checking whether filename already exists
            String patternStr1 = "^" + fileName.Replace(System.IO.Path.DirectorySeparatorChar, '%').Replace('/', '%') + @"_*[0-9]*\." + ext + @"$";
            String patternStr2 = fileName.Replace(System.IO.Path.DirectorySeparatorChar, '%').Replace('/', '%') + @"_*[0-9]*\." + ext + @"$";
            String temp;
            String dir = System.IO.Path.GetDirectoryName(fileName);
            dir = (dir == "" ? System.IO.Directory.GetCurrentDirectory() : dir);
            List<String> currFilesList = System.IO.Directory.GetFiles(dir, "*" + ext, System.IO.SearchOption.TopDirectoryOnly).ToList<String>();//.Select(e => e.Replace(System.IO.Path.DirectorySeparatorChar, '%')).ToList<String>();
            Int32 i = 0;
            while (i < currFilesList.Count)
            {
                currFilesList[i] = currFilesList[i].Replace(System.IO.Path.DirectorySeparatorChar, '%').Replace('/', '%');
                currFilesList[i] = System.Text.RegularExpressions.Regex.Match(currFilesList[i], patternStr2, System.Text.RegularExpressions.RegexOptions.Compiled).Value;
                i++;
            }
            try
            {
                i = 0;
                while (true)
                {
                    temp = currFilesList.First<String>(s => System.Text.RegularExpressions.Regex.IsMatch(s, patternStr1, System.Text.RegularExpressions.RegexOptions.Compiled));
                    if (!String.IsNullOrEmpty(temp))
                    {
                        currFilesList.Remove(temp);
                        i++;
                        continue;
                    }
                    else
                    {
                        break;
                    }
                }
            }
            catch (InvalidOperationException)
            {
                temp = (i == 0 ? String.Empty : "_" + i.ToString());
                fileName = fileName + temp + "." + ext;
            }
            return fileName;
        }

        public static void WriteColNamesTxtFile(System.IO.StreamWriter sw, String[] colNames)
        {
            String varColNames = "# " + colNames[0];
            for (int i = 1; i < colNames.Length; i++)
                varColNames += "\t" + colNames[i];
            sw.WriteLine(varColNames);
        }
        public static System.IO.StreamWriter CreateOutTxtFile(String fileName, String header)
        {
            fileName = KTzHeader.CheckAndGetFileName(fileName);
            System.IO.FileStream fs = new System.IO.FileStream(fileName, System.IO.FileMode.Create, System.IO.FileAccess.Write);
            System.IO.StreamWriter sw = new System.IO.StreamWriter(fs);
            Console.WriteLine("*** Writing to file: {0}", fileName);
            if (header[0] != '#')
                header = "#" + header;
            sw.WriteLine(header);
            sw.WriteLine("# {0}", DateTime.Now);
            sw.WriteLine("#-");
            return sw;
        }
        public static System.IO.StreamWriter CreateOutTxtFile(String fileName, String header, String[] colNames)
        {
            System.IO.StreamWriter sw = CreateOutTxtFile(fileName, header);
            WriteColNamesTxtFile(sw, colNames);
            return sw;
        }
        /// <summary>
        /// writes the columns to the supplied streamwriter
        /// </summary>
        /// <param name="sw">streamwriter object</param>
        /// <param name="colNames">name of the columns</param>
        /// <param name="c_data">c_data[i]: column i+1 in the file to be written</param>
        public static void WriteTxtFileColumns(System.IO.StreamWriter sw, String[] colNames, params Double[][] c_data)
        {
            if (colNames.Length == c_data.Length)
                colNames = (new String[] { "time" }).Concat(colNames).ToArray();
            WriteColNamesTxtFile(sw, colNames);
            WriteTxtFileColumns(sw, c_data);
        }

        public static void WriteTxtFileColumns(System.IO.StreamWriter sw, params Double[][] c_data)
        {
            int i;
            String outputStr = "{0}";
            for (i = 0; i < c_data.Length; i++)
                outputStr += "\t{" + (i + 1).ToString() + ":0.00000000e+000}";
            i = 0;
            //c_data = c_data.Transpose();
            while (i < c_data[0].Length)
            {
                sw.WriteLine(String.Format(outputStr, (new object[] { i }).Concat(c_data.GetCol(i).Cast<object>()).ToArray<object>()));
                i++;
            }
        }

        public static void CloseOutTxtFile(System.IO.StreamWriter sw)
        {
            sw.Close();
            Console.WriteLine("*** File written...");
            Console.WriteLine("-");
        }

        /// <summary>
        /// writes an output file of the required type (.dat for txt or .mat for mat)
        /// </summary>
        /// <typeparam name="TT">type of the data for the 1st and 2nd columns</typeparam>
        /// <typeparam name="T3">type of the data for the 3rd column</typeparam>
        /// <typeparam name="T4">type of the data for the 4th column</typeparam>
        /// <param name="param1">x-like param</param>
        /// <param name="param2">y-like param</param>
        /// <param name="z1">z-like values (it's first index must corresponde to the param1 index and it's second index must correspond to param2 index)</param>
        /// <param name="z1StdDev">standard deviation for z</param>
        /// <param name="param1Name">name of param1</param>
        /// <param name="param2Name">name of param2</param>
        /// <param name="z1Name">name of z</param>
        /// <param name="z1StdDevName">name of zStdDev</param>
        /// <param name="fileName">name of the file to be written</param>
        /// <param name="header">header of the file to be written</param>
        public static void Write11ColFile<TT, TU>(ref TT[] param1, ref TT[] param2,
            ref TU[][] z0, ref TU[][] z1, ref TU[][] z1StdDev, ref TU[][] z2, ref TU[][] z2StdDev, ref TU[][] z3, ref TU[][] z3StdDev, ref TU[][] z4, ref TU[][] z4StdDev,
            String[] colNames,
            String fileName, String header)
            where TT : struct
            where TU : struct
        {
            switch ((OutputFileFormat)KTzHeader.GetPar_Int32(KTzParameters.oFileFormat))
            {
                case OutputFileFormat.txt:
                    Write11ColFileTxt<TT, TU>(ref param1, ref param2, ref z0, ref z1, ref z1StdDev, ref z2, ref z2StdDev, ref z3, ref z3StdDev, ref z4, ref z4StdDev, colNames, fileName, header);
                    break;
                case OutputFileFormat.mat:
                    Write11ColFileMat<TT, TU>(ref param1, ref param2, ref z0, ref z1, ref z1StdDev, ref z2, ref z2StdDev, ref z3, ref z3StdDev, ref z4, ref z4StdDev, colNames, fileName, header);
                    break;
                default:
                    throw new ArgumentOutOfRangeException("oFileFormat", "unknown file format for output");
            }
        }

        /// <summary>
        /// writes an output file with 5 columns with the defined parameters
        /// </summary>
        /// <typeparam name="TT">type of the data for the 1st and 2nd columns</typeparam>
        /// <typeparam name="T3">type of the data for the 3rd column</typeparam>
        /// <typeparam name="T4">type of the data for the 4th column</typeparam>
        /// <param name="param1">x-like param</param>
        /// <param name="param2">y-like param</param>
        /// <param name="z1">z-like values (it's first index must corresponde to the param1 index and it's second index must correspond to param2 index)</param>
        /// <param name="z1StdDev">standard deviation for z</param>
        /// <param name="param1Name">name of param1</param>
        /// <param name="param2Name">name of param2</param>
        /// <param name="z1Name">name of z</param>
        /// <param name="z1StdDevName">name of zStdDev</param>
        /// <param name="fileName">name of the file to be written</param>
        /// <param name="ext">extension of the file to be written</param>
        /// <param name="header">header of the file to be written</param>
        public static void Write11ColFileTxt<TT, TU>(ref TT[] param1, ref TT[] param2,
            ref TU[][] z0, ref TU[][] z1, ref TU[][] z1StdDev, ref TU[][] z2, ref TU[][] z2StdDev, ref TU[][] z3, ref TU[][] z3StdDev, ref TU[][] z4, ref TU[][] z4StdDev,
            String[] colNames,
            String fileName, String header)
            where TT : struct
            where TU : struct

        {
            String ext = "dat";
            Int32 i, j, m, n;
            //System.IO.FileStream myFile;
            //System.IO.StreamWriter sw;

            // checking whether filename already exists
            fileName = KTzHeader.CheckAndGetFileName(fileName, ext);

            //Console.WriteLine("*** Writing to file: {0}", fileName);
            //myFile = new System.IO.FileStream(fileName, System.IO.FileMode.Create, System.IO.FileAccess.Write);
            //sw = new System.IO.StreamWriter(myFile);
            //sw.WriteLine(header);
            //sw.WriteLine("# {0}", DateTime.Now);
            //sw.WriteLine("#-");

            String fmtP = (typeof(TT) == typeof(Double) ? "0.00000000e+000" : String.Empty);
            String fmtV = (typeof(TU) == typeof(Double) ? "0.00000000e+000" : String.Empty);
            //String outputStr = "{0:" + fmtP + "}\t{1:" + fmtP + "}\t{2:" + fmtV + "}\t{3:" + fmtV + "}\t{4:" + fmtV + "}\t{5:" + fmtV + "}\t{6:" + fmtV + "}\t{7:" + fmtV + "}\t{8:" + fmtV + "}";
            //String varColNames = colNames[2] + "\t" + colNames[3] + "\t" + colNames[4] + "\t" + colNames[5] + "\t" + colNames[6] + "\t" + colNames[7] + "\t" + colNames[8];
            String outputStr = "{0:" + fmtP + "}\t{1:" + fmtP + "}\t{2:" + fmtV + "}";
            for (i = 3; i < 11; i++)
            {
                outputStr += "\t{" + i.ToString() + ":" + fmtV + "}";
            }

            // setting the bounds of the loops
            // this condition is necessary because I'll never know whether it is param1 or param2 which correspond to the first dimension of z1 and z2
            TT[] p1, p2;
            if (z1.Length == param2.Length)
            {
                m = param2.Length;
                n = param1.Length;
                p1 = param2;
                p2 = param1;
                //sw.WriteLine("#" + colNames[1] + "\t" + colNames[0] + "\t" + varColNames);
            }
            else
            {
                m = param1.Length;
                n = param2.Length;
                p1 = param1;
                p2 = param2;
                String temp = colNames[1];
                colNames[1] = colNames[0];
                colNames[0] = temp;
                //sw.WriteLine("#" + colNames[0] + "\t" + colNames[1] + "\t" + varColNames);
            }

            System.IO.StreamWriter sw = KTzHeader.CreateOutTxtFile(fileName, header, colNames);

            // writing data
            i = 0;
            while (i < m)
            {
                j = 0;
                while (j < n)
                {
                    sw.WriteLine(outputStr, p1[i], p2[j], z0[i][j], z1[i][j], z1StdDev[i][j], z2[i][j], z2StdDev[i][j], z3[i][j], z3StdDev[i][j], z4[i][j], z4StdDev[i][j]);
                    j++;
                }
                i++;
            }


            // closing file
            //sw.Close();
            //myFile.Close();
            KTzHeader.CloseOutTxtFile(sw);
        }


        public static void Write11ColFileMat<TT, TU>(ref TT[] param1, ref TT[] param2,
            ref TU[][] z0, ref TU[][] z1, ref TU[][] z1StdDev, ref TU[][] z2, ref TU[][] z2StdDev, ref TU[][] z3, ref TU[][] z3StdDev, ref TU[][] z4, ref TU[][] z4StdDev,
            String[] colNames, String fileName, String header)
            where TT : struct
            where TU : struct
        {
            String ext = "mat";

            // checking whether filename already exists
            fileName = KTzHeader.CheckAndGetFileName(fileName, ext);
            Console.WriteLine("*** Writing to file: {0}", fileName);

            // setting the bounds of the loops
            // this condition is necessary because I'll never know whether it is param1 or param2 which correspond to the first dimension of z1 and z2
            TT[] p1, p2;
            if (z1.Length == param2.Length)
            {
                p1 = param2;
                p2 = param1;
                //sw.WriteLine("#" + colNames[1] + "\t" + colNames[0] + "\t" + varColNames);
            }
            else
            {
                p1 = param1;
                p2 = param2;
                String temp = colNames[1];
                colNames[1] = colNames[0];
                colNames[0] = temp;
                //sw.WriteLine("#" + colNames[0] + "\t" + colNames[1] + "\t" + varColNames);
            }


            colNames = KTzHeader.KeepLetterOrDigit(colNames);

            var matDataBuilder = new MatFileHandler.DataBuilder();
            Dictionary<KTzParamGroups, MatFileHandler.IVariable> simParams = KTzHeader.GetAllParamPairsAsMatlabStruct(matDataBuilder);

            var p1_mat = matDataBuilder.NewArray(p1, p1.Length, 1);
            var p2_mat = matDataBuilder.NewArray(p2, 1, p2.Length);
            var z0_mat = matDataBuilder.NewArray(z0, z0.Length, z0[0].Length);
            var z1_mat = matDataBuilder.NewArray(z1, z1.Length, z1[0].Length);
            var z2_mat = matDataBuilder.NewArray(z2, z2.Length, z2[0].Length);
            var z3_mat = matDataBuilder.NewArray(z3, z3.Length, z3[0].Length);
            var z4_mat = matDataBuilder.NewArray(z4, z4.Length, z4[0].Length);
            var z1StdDev_mat = matDataBuilder.NewArray(z1StdDev, z1StdDev.Length, z1StdDev[0].Length);
            var z2StdDev_mat = matDataBuilder.NewArray(z2StdDev, z2StdDev.Length, z2StdDev[0].Length);
            var z3StdDev_mat = matDataBuilder.NewArray(z3StdDev, z3StdDev.Length, z3StdDev[0].Length);
            var z4StdDev_mat = matDataBuilder.NewArray(z4StdDev, z4StdDev.Length, z4StdDev[0].Length);
            using (var fileStream = new System.IO.FileStream(fileName, System.IO.FileMode.Create))
            {
                var writer = new MatFileHandler.MatFileWriter(fileStream);
                writer.Write(matDataBuilder.NewFile(simParams.Values.Concat(new[]
                { matDataBuilder.NewVariable(colNames[0] , p1_mat),
                  matDataBuilder.NewVariable(colNames[1] , p2_mat),
                  matDataBuilder.NewVariable(colNames[2] , z0_mat),
                  matDataBuilder.NewVariable(colNames[3] , z1_mat),
                  matDataBuilder.NewVariable(colNames[4] , z1StdDev_mat),
                  matDataBuilder.NewVariable(colNames[5] , z2_mat),
                  matDataBuilder.NewVariable(colNames[6] , z2StdDev_mat),
                  matDataBuilder.NewVariable(colNames[7] , z3_mat),
                  matDataBuilder.NewVariable(colNames[8] , z3StdDev_mat),
                  matDataBuilder.NewVariable(colNames[9] , z4_mat),
                  matDataBuilder.NewVariable(colNames[10], z4StdDev_mat),
                  matDataBuilder.NewVariable("file_header", matDataBuilder.NewCharArray(header)),
                  matDataBuilder.NewVariable("file_info", matDataBuilder.NewCharArray(String.Format("data matrices: rows -> {0}; cols -> {1}",colNames[0],colNames[1])))
                })));
            }
            /**/
            /*
            System.IO.StreamWriter sw = KTzHeader.CreateOFile(ref fileName, header, colNames);

            // writing data
            i = 0;
            while (i < m)
            {
                j = 0;
                while (j < n)
                {
                    sw.WriteLine(outputStr, p1[i], p2[j], z0[i][j], z1[i][j], z1StdDev[i][j], z2[i][j], z2StdDev[i][j], z3[i][j], z3StdDev[i][j], z4[i][j], z4StdDev[i][j]);
                    j++;
                }
                i++;
            }

            // closing file
            //sw.Close();
            //myFile.Close();
            KTzHeader.CloseOFile(sw);
            /**/
        }

        public static String[] KeepLetterOrDigit(String[] str)
        {
            for (var i = 0; i < str.Length; i++)
            {
                str[i] = KeepLetterOrDigit(str[i]);
            }
            return str;
        }

        public static String KeepLetterOrDigit(String str)
        {
            String s = new String(str.Where(c => char.IsLetterOrDigit(c) || c == '_').ToArray());
            var rgx = new System.Text.RegularExpressions.Regex("^[\\d_]+");
            return rgx.Replace(s,"");
        }

        /// <summary>
        /// get the executable name - useful for the usage of this program through command line
        /// </summary>
        /// <returns>the executable name of the current program</returns>
        public static String GetKTzV2ExeName()
        {
            return System.IO.Path.GetFileName(Environment.GetCommandLineArgs()[0]);
        }
    }

    public enum KTzParameters
    {
        /// <summary>
        /// neuron parameter K
        /// </summary>
        K,
        /// <summary>
        /// neuron parameter T
        /// </summary>
        T,
        /// <summary>
        /// neuron parameter delta - refractory period and burst damping
        /// </summary>
        d,
        /// <summary>
        /// neuron parameter lambda - refractory period
        /// </summary>
        l,
        /// <summary>
        /// neuron parameter x_R - burst duration
        /// </summary>
        xR,
        /// <summary>
        /// neuron parameter H - external field, polarizing current
        /// </summary>
        H,
        /// <summary>
        /// neuron parameter Q KTz2Tanh X variable - external field, polarizing current
        /// </summary>
        Q,
        /// <summary>
        /// neuron initial condition on x
        /// </summary>
        x0,
        /// <summary>
        /// neuron initial condition on y
        /// </summary>
        y0,
        /// <summary>
        /// neuron initial condition on z
        /// </summary>
        z0,
        /// <summary>
        /// synapse parameter tau_f - time delay on f variable
        /// </summary>
        tauf,
        /// <summary>
        /// synapse parameter tau_g - time delay on g variable
        /// </summary>
        taug,
        /// <summary>
        /// synapse parameter - coupling intensity
        /// </summary>
        J,
        /// <summary>
        /// synapse parameter - minimum J on bifurcation simulations
        /// </summary>
        minJ,
        /// <summary>
        /// synapse parameter - maximum J on bifurcation simulations
        /// </summary>
        maxJ,
        /// <summary>
        /// synapse parameter - range values for J parameter when running Bifurcation simulation (uses MATLAB syntax)
        /// </summary>
        JRange,
        /// <summary>
        /// stimulus parameter - stimulus intensity
        /// </summary>
        I,
        /// <summary>
        /// stimulus parameter - minimum I on bifurcation simulations
        /// </summary>
        minI,
        /// <summary>
        /// stimulus parameter - maximum I on bifurcation simulations
        /// </summary>
        maxI,
        /// <summary>
        /// synapse parameter - range values for I parameter when running Bifurcation simulation (uses MATLAB syntax)
        /// </summary>
        IRange,
        /// <summary>
        /// stimulus parameter - std deviation of I in case of Poisson stimulus
        /// </summary>
        IStdDev,
        /// <summary>
        /// Poisson stimulus parameter - stimulus rate
        /// </summary>
        r,
        /// <summary>
        /// Poisson stimulus parameter - minimum stimulus rate on bifurcation simulations
        /// </summary>
        minr,
        /// <summary>
        /// Poisson stimulus parameter - maximum stimulus rate on bifurcation simulations
        /// </summary>
        maxr,
        /// <summary>
        /// synapse parameter - range values for r parameter when running Bifurcation simulation (uses MATLAB syntax)
        /// </summary>
        rRange,
        /// <summary>
        /// Noisy synapse parameter - noise amplitude
        /// </summary>
        R,
        /// <summary>
        /// Noisy synapse parameter - ratio of noise in case of ProportionalAmplitude noisy synapse
        /// </summary>
        noiseRatio,
        /// <summary>
        /// Dynamical synapse parameter - coupling intensity * u
        /// </summary>
        alpha,
        /// <summary>
        /// Dynamical synapse parameter - minimum coupling intensity on bifurcation simulations
        /// </summary>
        minalpha,
        /// <summary>
        /// Dynamical synapse parameter - maximum coupling intensity on bifurcation simulations
        /// </summary>
        maxalpha,
        /// <summary>
        /// synapse parameter - range values for alpha parameter when running Bifurcation simulation (uses MATLAB syntax)
        /// </summary>
        alphaRange,
        /// <summary>
        /// Dynamical synapse parameter - decay of the synapse
        /// </summary>
        u,
        /// <summary>
        /// Dynamical synapse parameter - delay on dynamical coupling J
        /// </summary>
        tauJ,
        /// <summary>
        /// Dynamical synapse parameter - precision on integrating J equations (with Euler method)
        /// </summary>
        dt,
        /// <summary>
        /// Watts-Strogatz or Random Graph Network parameter - rewire probability or edge creation prob
        /// </summary>
        rewP,
        /// <summary>
        /// Sampling paramemter - fraction of neurons to sample on the network
        /// </summary>
        sampFrac,
        /// <summary>
        /// Network parameter - Amount of elements on x axis (linear size) (total amount is always assumed as Lx*Ly*Lz)
        /// </summary>
        Lx,
        /// <summary>
        /// Network parameter - Amount of elements on y axis (linear size) (total amount is always assumed as Lx*Ly*Lz)
        /// </summary>
        Ly,
        /// <summary>
        /// Network parameter - Amount of elements on z axis (linear size) (total amount is always assumed as Lx*Ly*Lz)
        /// </summary>
        Lz,
        /// <summary>
        /// Simulation parameter - total timesteps to run each simulation
        /// </summary>
        nSteps,
        /// <summary>
        /// Simulation parameter - if "ProportionalToPoissonRate" then (nSteps - nStartStep) = 10/(r*N), where r-> Poisson rate; N=Lx*Ly*Lz (number of neurons)
        /// </summary>
        simTimeScheme,
        /// <summary>
        /// Simulation parameter - timestep in which the network starts to be recorded
        /// </summary>
        nStart,
        /// <summary>
        /// Network parameter - type of the network
        /// </summary>
        netType,
        /// <summary>
        /// Network parameter - directed = yes: the synapses are one-way i->j; directed = no: two-way synapses
        /// </summary>
        netDir,
        /// <summary>
        /// Synapse parameter - type of the synapse
        /// </summary>
        sType,
        /// <summary>
        /// Synapse parameter - synaptic noise type
        /// </summary>
        noiseType,
        /// <summary>
        /// Simulation parameter - type of initial condition
        /// </summary>
        iCond,
        /// <summary>
        /// Network parameter - dimension of the network (network will have L^dim neurons)
        /// </summary>
        dim,
        /// <summary>
        /// Stimulus parameter - type of the stimulus
        /// </summary>
        stimType,
        /// <summary>
        /// Stimulus parameter - timestep in which the network will be stimulated
        /// </summary>
        sStim,
        /// <summary>
        /// Stimulus parameter - index of the neuron which will be stimulated
        /// </summary>
        iStim,
        /// <summary>
        /// stimulus paremeter - method for choosing stimulus index: Fixed, Random or MostConnected neuron
        /// </summary>
        indChoice,
        /// <summary>
        /// Simulation parameter - size of the bin (in timesteps) to divide nSteps in order to run the simulation
        /// </summary>
        tBin,
        /// <summary>
        /// Simulation parameter - choose whether the program will run a simulation only to write output data in the case of dynamics simulation
        /// </summary>
        wData,
        /// <summary>
        /// if wData=Yes, then wCSV=Yes writes a CSV file with the x_i(t) data
        /// </summary>
        wCSV,
        /// <summary>
        /// Simulation parameter - choose whether the program will run a simulation only to write an avalanche spike distribution file in the case of dynamics simulation
        /// </summary>
        wAvalDist,
        /// <summary>
        /// Simulation parameter - chooses whether to save the spike times for each neuron in a simulation
        /// </summary>
        saveSpikeTimes,
        /// <summary>
        /// Debug parameter - write observations of the avalanche sizes for Bifurcation sim type for Debug purposes
        /// </summary>
        wObs,
        /// <summary>
        /// Simulation parameter - choose whether the program will write a file containing xi-xj data
        /// </summary>
        wDif,
        /// <summary>
        /// Simulation parameter - amount of J on bifurcation simulations
        /// </summary>
        nJ,
        /// <summary>
        /// Simulation parameter - amount of I on bifurcation simulations
        /// </summary>
        nI,
        /// <summary>
        /// Simulation parameter - amount of r on bifurcation simulations
        /// </summary>
        nr,
        /// <summary>
        /// Simulation parameter - amount of alpha on bifurcation simulations
        /// </summary>
        nalpha,
        /// <summary>
        /// Simulation parameter - amount of realizations for each pair (stimulus,coupling) on bifurcation simulations
        /// </summary>
        nSim,
        /// <summary>
        /// Delta Stimulus parameter - time interval between two consecutive delta stimuli
        /// </summary>
        deltaT,
        /// <summary>
        /// Watts-Strogatz network parameter - amount of neighbours in the initial configuration of WS network
        /// </summary>
        nNeigh,
        /// <summary>
        /// Barabasi-Albert network parameter - amount of connections which a new node will initially have when creating BA network
        /// </summary>
        nConn,
        /// <summary>
        /// Simulation parameter - choose if the program will count number of neurons spiking or number of spikes
        /// </summary>
        cVar,
        /// <summary>
        /// Neuron parameter - choose whether the neuron will average its synaptic input (may be useful in mean-field networks)
        /// </summary>
        avgInp,
        /// <summary>
        /// Simulation parameter - choose the sampling type (full or subsampled)
        /// </summary>
        samp,
        /// <summary>
        /// No activity Stimulus parameter - interval (in tBin's) between an avalanche and a new stimulus for no activity stimulus type
        /// </summary>
        rest,
        /// <summary>
        /// Simulation parameter - choose if the program will run bifurcation or dynamics
        /// </summary>
        simType,
        /// <summary>
        /// Simulation parameter - each avalanche is generated independently of each other or after a time window
        /// </summary>
        dynType,
        /// <summary>
        /// Simulation parameter - prefix of the output file name
        /// </summary>
        oFile,
        /// <summary>
        /// chooses whether to write to txt or mat files
        /// </summary>
        oFileFormat,
        /// <summary>
        /// Simulation parameter - input file with adjacency matrix
        /// </summary>
        netFile,
        /// <summary>
        /// Neuron parameter - type of neuron to use
        /// </summary>
        neuron,
        /// <summary>
        /// simulation parameter - OnTheFly == writes data after each pair (I,J); OnTheEnd == writes after everything is finished
        /// </summary>
        bifWrite,
        /// <summary>
        /// simulation parameter - determines whether the output average should be taken over many avalanches or just temporal average (for bifurcation simulations)
        /// </summary>
        outAvgMode,
        /// <summary>
        /// simulation parameter - writes rho(t) for each J in JRange (only if outAvgMode=OverTime); rho(t) = sum x(t)>0
        /// </summary>
        writeRhoTS,
        /// <summary>
        /// neuron parameter - threshold of the SI(RS) element
        /// </summary>
        Theta,
        /// <summary>
        /// coupling parameter - determines whether the J are given by J parameter or by adjacency matrix
        /// </summary>
        coupParam,
        /// <summary>
        /// sim parameter - min value for threshold of SI(RS) element in bifurcation simulation
        /// </summary>
        minTheta,
        /// <summary>
        /// sim parameter - max value for threshold of SI(RS) element in bifurcation simulation
        /// </summary>
        maxTheta,
        /// <summary>
        /// sim parameter - number of values for threshold of SI(RS) element in bifurcation simulation
        /// </summary>
        nTheta,
        /// <summary>
        /// sim parameter - range of values for threshold of SI(RS) element in bifurcation simulation (uses MATLAB syntax)
        /// </summary>
        ThetaRange,
        /// <summary>
        /// disorder for K
        /// </summary>
        KDisorder,
        /// <summary>
        /// disorder for T
        /// </summary>
        TDisorder,
        /// <summary>
        /// disorder for d
        /// </summary>
        dDisorder,
        /// <summary>
        /// disorder for l
        /// </summary>
        lDisorder,
        /// <summary>
        /// disorder for xR
        /// </summary>
        xRDisorder,
        /// <summary>
        /// disorder for H
        /// </summary>
        HDisorder,
        /// <summary>
        /// disorder for Q
        /// </summary>
        QDisorder,
        /// <summary>
        /// disorder for J
        /// </summary>
        JDisorder,
        /// <summary>
        /// disorder for tauf
        /// </summary>
        taufDisorder,
        /// <summary>
        /// disorder for taug
        /// </summary>
        taugDisorder,
        /// <summary>
        /// disorder for R
        /// </summary>
        RDisorder,
        /// <summary>
        /// disorder for noiseRatio
        /// </summary>
        noiseRatioDisorder,
        /// <summary>
        /// disorder for alpha
        /// </summary>
        alphaDisorder,
        /// <summary>
        /// disorder for u
        /// </summary>
        uDisorder,
        /// <summary>
        /// disorder for tauJ
        /// </summary>
        tauJDisorder,
        /// <summary>
        /// disorder for Theta
        /// </summary>
        ThetaDisorder,
        /// <summary>
        /// parameter used to simulate a range in DynamicsWithinParamRange simulation
        /// </summary>
        ParamForRange
    }
}
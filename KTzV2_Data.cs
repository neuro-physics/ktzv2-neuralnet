using System;
using System.Collections.Generic;
using System.IO;
using System.Text.RegularExpressions;
using System.Linq;
using System.Xml;
using KTzV2.Maths.Random;
using MatFileHandler;

namespace KTzV2.Data
{
    public class ThermoStatistics
    {
        public Double mMax { get; private set; }
        public Double m { get; private set; }
        public Double mStdDev { get; private set; }
        public Double chi { get; private set; }
        public Double chiStdDev { get; private set; }
        public Double U4 { get; private set; }
        public Double U4StdDev { get; private set; }
        public Double[] mi { get; private set; }
        /// <summary>
        /// calculates order param, susceptibility and 4th order Binder cumulants, given  a set mi of measurments of the order parameter; variances are estimated using jackknife
        /// </summary>
        /// <param name="m">average order param</param>
        /// <param name="m2">average of squared order param</param>
        /// <param name="m4">average of 4th power of order param</param>
        /// <param name="mi">set of measurements of order param</param>
        /// <param name="miMax">max value wihtin mi set</param>
        /// <param name="N">amount of neurons or spins in the network</param>
        public ThermoStatistics(Double m, Double m2, Double m4, Double[] mi, Double miMax, Int32 N)
        {
            //m[i][j] /= nSim; // calculating mean spikeCountTot - as the mean variable is an Int32, I'm losing 1 spike doing this division
            //chi[i][j] = nNeurons * (mStdDev[i][j] / nSim - m[i][j] * m[i][j]);
            //U4[i][j] = 1 - (U4[i][j] / nSim) / (3 * (mStdDev[i][j] / nSim) * (mStdDev[i][j] / nSim));
            //mStdDev[i][j] = Math.Sqrt(mStdDev[i][j] / nSim - m[i][j] * m[i][j]); // mean of square minus square of mean
            this.mMax = miMax;
            this.m = m;
            this.mStdDev = Math.Sqrt(m2 - m * m);
            this.chi = N * (m2 - m * m);
            this.U4 = 1.0D - m4 / (3.0D * m2 * m2);
            this.mi = mi;
            this.JackKnife(mi);
        }

        public void WriteObservations(String fileName, String header, KTzV2.Data.Header.OutputFileFormat outFileFormat = KTzV2.Data.Header.OutputFileFormat.mat)
        {
            switch (outFileFormat)
            {
                case KTzV2.Data.Header.OutputFileFormat.txt:
                    this.WriteObservationsTxt(fileName, header);
                    break;
                case KTzV2.Data.Header.OutputFileFormat.mat:
                    this.WriteObservationsMat(fileName, header);
                    break;
                default:
                    throw new ArgumentOutOfRangeException("outFileFormat", "unknown file format for output");
            }
        }

        private void WriteObservationsTxt(String fileName, String header)
        {
            FileStream fs = new FileStream(fileName, FileMode.Create, FileAccess.Write);
            StreamWriter sw = new StreamWriter(fs);
            sw.WriteLine(header);
            for (int i = 0; i < this.mi.Length; i++)
            {
                sw.WriteLine(this.mi[i]);
            }
            sw.Close();
            fs.Close();
        }
        private void WriteObservationsMat(String fileName, String header)
        {
            using (var fileStream = new System.IO.FileStream(fileName, System.IO.FileMode.Create))
            {
                var matDataBuilder = new MatFileHandler.DataBuilder();
                var writer = new MatFileHandler.MatFileWriter(fileStream);
                writer.Write(
                    matDataBuilder.NewFile(new[]{
                         matDataBuilder.NewVariable("order_param_observations", matDataBuilder.NewArray(this.mi, this.mi.Length, 1)),
                         matDataBuilder.NewVariable("file_header", matDataBuilder.NewCharArray(header))
                    })
                );
            }
        }
        private void JackKnife(Double[] mi)
        {
            int i = 0, j;
            Double c, c2, u, u2, ss, ss2, ss4, cM, uM;
            int N = mi.Length;
            c = 0.0D;
            c2 = 0.0D;
            u = 0.0D;
            u2 = 0.0D;
            while (i < N)
            {
                ss = 0.0D;
                ss2 = 0.0D;
                ss4 = 0.0D;
                j = 0;
                while (j < i)
                {
                    ss += mi[j];
                    ss2 += mi[j] * mi[j];
                    ss4 += mi[j] * mi[j] * mi[j] * mi[j];
                    j++;
                }
                j++;
                while (j < N)
                {
                    ss += mi[j];
                    ss2 += mi[j] * mi[j];
                    ss4 += mi[j] * mi[j] * mi[j] * mi[j];
                    j++;
                }
                ss = ss / (Double)(N - 1);
                ss2 = ss2 / (Double)(N - 1);
                cM = ss2 - ss * ss;
                c += cM;
                c2 += cM * cM;
                uM = 1.0D - ((Double)ss4 / (Double)(N - 1)) / (3.0D * ss2 * ss2);
                u += uM;
                u2 += uM * uM;
                i++;
            }
            c = c / (Double)N;
            u = u / (Double)N;
            this.chiStdDev = Math.Sqrt(c2 / (Double)N - c * c);
            this.U4StdDev = Math.Sqrt(u2 / (Double)N - u * u);
        }
    }

    public enum NetworkSamplingType
    {
        Full,
        Partial
    }

    public enum CountVariable
    {
        NumberOfNeurons,
        NumberOfSpikes
    }

    public enum YesOrNoAnswer
    {
        No,
        Yes
    }

    public enum InitialConditionType
    {
        FromXMLFile,
        FromLastSimulation,
        ProgramSpecified,
        RandomXPositively,
        RandomYPositively,
        RandomZPositively,
        RandomXYPositively,
        RandomXZPositively,
        RandomYZPositively,
        RandomXYZPositively,
        RandomX,
        RandomY,
        RandomZ,
        RandomXY,
        RandomXZ,
        RandomYZ,
        RandomXYZ
    }

    public class FGetInitCond
    {
        public static GetICInstance use(InitialConditionType icType, Int32 nSims, Int32 nNeurons)
        {
            switch (icType)
            {
                case InitialConditionType.FromXMLFile: return new GetICFromXMLFile(nSims, nNeurons);
                case InitialConditionType.FromLastSimulation: return new GetICLastSim();
                case InitialConditionType.ProgramSpecified: return new GetICSpec();
                case InitialConditionType.RandomXPositively: return new GetICRandXPos();
                case InitialConditionType.RandomYPositively: return new GetICRandYPos();
                case InitialConditionType.RandomZPositively: return new GetICRandZPos();
                case InitialConditionType.RandomXYPositively: return new GetICRandXYPos();
                case InitialConditionType.RandomXZPositively: return new GetICRandXZPos();
                case InitialConditionType.RandomYZPositively: return new GetICRandYZPos();
                case InitialConditionType.RandomXYZPositively: return new GetICRandXYZPos();
                case InitialConditionType.RandomX: return new GetICRandX();
                case InitialConditionType.RandomY: return new GetICRandY();
                case InitialConditionType.RandomZ: return new GetICRandZ();
                case InitialConditionType.RandomXY: return new GetICRandXY();
                case InitialConditionType.RandomXZ: return new GetICRandXZ();
                case InitialConditionType.RandomYZ: return new GetICRandYZ();
                case InitialConditionType.RandomXYZ: return new GetICRandXYZ();
            }
            throw new ArgumentOutOfRangeException("Initial condition unknown!");
        }
    }

    public abstract class IGetInitCond
    {
        //protected HomogeneousRand rand { get; set; }
        public virtual Double getX(Double xi) { throw new NotImplementedException(); }
        public virtual Double getY(Double yi) { throw new NotImplementedException(); }
        public virtual Double getZ(Double zi) { throw new NotImplementedException(); }
        public virtual Double getX(Int32 n)
        {
            Double temp = 0.0;
            while (true)
            {
                Console.Write("x{0} = ", n);
                try
                {
                    temp = Convert.ToDouble(Console.ReadLine());
                    break;
                }
                catch (Exception e)
                {
                    Console.WriteLine(e.Message);
                    continue;
                }
            }
            return temp;
        }
        public virtual Double getY(Int32 n)
        {
            Double temp = 0.0;
            while (true)
            {
                Console.Write("y{0} = ", n);
                try
                {
                    temp = Convert.ToDouble(Console.ReadLine());
                    break;
                }
                catch (Exception e)
                {
                    Console.WriteLine(e.Message);
                    continue;
                }
            }
            return temp;
        }
        public virtual Double getZ(Int32 n)
        {
            Double temp = 0.0;
            while (true)
            {
                Console.Write("z{0} = ", n);
                try
                {
                    temp = Convert.ToDouble(Console.ReadLine());
                    break;
                }
                catch (Exception e)
                {
                    Console.WriteLine(e.Message);
                    continue;
                }
            }
            return temp;
        }
        public virtual Double getX(Int32 nSim, Int32 nNeu) { throw new NotImplementedException(); }
        public virtual Double getY(Int32 nSim, Int32 nNeu) { throw new NotImplementedException(); }
        public virtual Double getZ(Int32 nSim, Int32 nNeu) { throw new NotImplementedException(); }
    }

    public class GetICInstance : IGetInitCond
    {
        public GetICInstance() { }
    }

    public class GetICFromXMLFile : GetICInstance
    {
        private Double[][] x0 { get; set; }
        private Double[][] y0 { get; set; }
        private Double[][] z0 { get; set; }
        public String fileName { get; private set; }

        public GetICFromXMLFile(Int32 nSims, Int32 nNeurons)
        {
            Console.WriteLine("*** Searching XML IC files...");
            String patternStr = "ic_ns" + nSims + "_nn" + nNeurons + @"_ty-*[0-9]*_*[0-9]*\.xml$";
            //String patternStr2 = "^ic_ns" + nSims + "_nn" + nNeurons + @"_ty-*[0-9]*_*[0-9]*\.xml$";
            String currDir = Directory.GetCurrentDirectory();
            List<String> availableFiles = Directory.GetFiles(currDir, "*.xml", SearchOption.TopDirectoryOnly).ToList<String>();
            availableFiles.RemoveAll(s => !Regex.IsMatch(s, patternStr, RegexOptions.Compiled));
            Int32 n = availableFiles.Count;
            if (n > 0)
            {
                Int32 i, k;
                if (n > 1)
                {
                    i = 0;
                    Console.WriteLine("*** Choose one of the following files to load:");
                    while (i < n)
                    {
                        availableFiles[i] = Regex.Match(availableFiles[i], patternStr, RegexOptions.Compiled).Value;
                        Console.WriteLine("*** {0}: {1}", i, availableFiles[i]);
                        i++;
                    }
                    k = -1;
                    while (k < 0 || k >= n)
                    {
                        Console.Write("*** file index = ");
                        try
                        {
                            k = Convert.ToInt32(Console.ReadLine());
                        }
                        catch (Exception e)
                        {
                            Console.WriteLine("*** index error: {0}", e.Message);
                            k = -1;
                            continue;
                        }
                    }
                }
                else
                {
                    k = 0;
                }
                Console.WriteLine("*** Loading XML IC file: {0}", availableFiles[k]);
                this.fileName = currDir + Path.DirectorySeparatorChar + availableFiles[k];
                XmlDocument xmlDoc = new XmlDocument();
                xmlDoc.Load(this.fileName);

                x0 = new Double[nSims][];
                y0 = new Double[nSims][];
                z0 = new Double[nSims][];

                i = 0;
                Int32 j;
                XmlNode simNode;
                XmlNodeList x0List, y0List, z0List;
                while (i < nSims)
                {
                    simNode = xmlDoc.SelectSingleNode("//Simulation[@order='" + (i + 1).ToString() + "']");
                    x0List = simNode.SelectNodes("//x0");
                    y0List = simNode.SelectNodes("//y0");
                    z0List = simNode.SelectNodes("//z0");
                    x0[i] = new Double[nNeurons];
                    y0[i] = new Double[nNeurons];
                    z0[i] = new Double[nNeurons];
                    j = 0;
                    while (j < nNeurons)
                    {
                        x0[i][j] = Convert.ToDouble(x0List[j].InnerText);
                        y0[i][j] = Convert.ToDouble(y0List[j].InnerText);
                        z0[i][j] = Convert.ToDouble(z0List[j].InnerText);
                        j++;
                    }
                    i++;
                }
                Console.WriteLine("*** File loaded successfully!");
            }
            else
            {
                throw new Exception("No file available!");
            }
        }
        public override Double getX(Int32 nSim, Int32 nNeu) { return x0[nSim][nNeu]; }
        public override Double getY(Int32 nSim, Int32 nNeu) { return y0[nSim][nNeu]; }
        public override Double getZ(Int32 nSim, Int32 nNeu) { return z0[nSim][nNeu]; }
    }

    public class GetICRandXPos : GetICInstance
    {
        private HomogeneousRand rand { get; set; }

        public GetICRandXPos() { rand = new HomogeneousRand(); }

        public override Double getX(Double xi)
        {
            return xi * rand.GetRandomFull();
        }
        public override Double getY(Double yi) { return yi; }
        public override Double getZ(Double zi) { return zi; }
    }

    public class GetICRandYPos : GetICInstance
    {
        private HomogeneousRand rand { get; set; }

        public GetICRandYPos() { rand = new HomogeneousRand(); }

        public override Double getY(Double yi)
        {
            return yi * rand.GetRandomFull();
        }
        public override Double getX(Double xi) { return xi; }
        public override Double getZ(Double zi) { return zi; }
    }

    public class GetICRandZPos : GetICInstance
    {
        private HomogeneousRand rand { get; set; }

        public GetICRandZPos() { rand = new HomogeneousRand(); }

        public override Double getZ(Double zi)
        {
            return zi * rand.GetRandomFull();
        }
        public override Double getY(Double yi) { return yi; }
        public override Double getX(Double xi) { return xi; }
    }

    public class GetICRandXYPos : GetICInstance
    {
        private HomogeneousRand rand { get; set; }

        public GetICRandXYPos() { rand = new HomogeneousRand(); }

        public override Double getX(Double xi)
        {
            return xi * rand.GetRandomFull();
        }
        public override Double getY(Double yi)
        {
            return yi * rand.GetRandomFull();
        }
        public override Double getZ(Double zi) { return zi; }
    }

    public class GetICRandXZPos : GetICInstance
    {
        private HomogeneousRand rand { get; set; }

        public GetICRandXZPos() { rand = new HomogeneousRand(); }

        public override Double getX(Double xi)
        {
            return xi * rand.GetRandomFull();
        }
        public override Double getZ(Double zi)
        {
            return zi * rand.GetRandomFull();
        }
        public override Double getY(Double yi) { return yi; }
    }

    public class GetICRandYZPos : GetICInstance
    {
        private HomogeneousRand rand { get; set; }

        public GetICRandYZPos() { rand = new HomogeneousRand(); }

        public override Double getY(Double yi)
        {
            return yi * rand.GetRandomFull();
        }
        public override Double getZ(Double zi)
        {
            return zi * rand.GetRandomFull();
        }
        public override Double getX(Double xi) { return xi; }
    }

    public class GetICRandXYZPos : GetICInstance
    {
        private HomogeneousRand rand { get; set; }

        public GetICRandXYZPos() { rand = new HomogeneousRand(); }

        public override Double getX(Double xi)
        {
            return xi * rand.GetRandomFull();
        }
        public override Double getY(Double yi)
        {
            return yi * rand.GetRandomFull();
        }
        public override Double getZ(Double zi)
        {
            return zi * rand.GetRandomFull();
        }
    }

    public class GetICLastSim : GetICInstance
    {
        private HomogeneousRand rand { get; set; }
        public GetICLastSim() { }
        public override Double getX(Double xi) { return xi; }
        public override Double getY(Double yi) { return yi; }
        public override Double getZ(Double zi) { return zi; }
    }

    public class GetICSpec : GetICInstance
    {
        private HomogeneousRand rand { get; set; }
        public GetICSpec() { }
        public override Double getX(Double xi) { return xi; }
        public override Double getY(Double yi) { return yi; }
        public override Double getZ(Double zi) { return zi; }
    }

    public class GetICRandX : GetICInstance
    {
        private HomogeneousRand rand { get; set; }

        public GetICRandX() { rand = new HomogeneousRand(); }

        public override Double getX(Double xi)
        {
            return xi * (rand.GetRandomFull() * 2.0 - 1.0);
        }
        public override Double getY(Double yi) { return yi; }
        public override Double getZ(Double zi) { return zi; }
    }

    public class GetICRandY : GetICInstance
    {
        private HomogeneousRand rand { get; set; }

        public GetICRandY() { rand = new HomogeneousRand(); }

        public override Double getY(Double yi)
        {
            return yi * (rand.GetRandomFull() * 2.0 - 1.0);
        }
        public override Double getX(Double xi) { return xi; }
        public override Double getZ(Double zi) { return zi; }
    }

    public class GetICRandZ : GetICInstance
    {
        private HomogeneousRand rand { get; set; }

        public GetICRandZ() { rand = new HomogeneousRand(); }

        public override Double getZ(Double zi)
        {
            return zi * (rand.GetRandomFull() * 2.0 - 1.0);
        }
        public override Double getY(Double yi) { return yi; }
        public override Double getX(Double xi) { return xi; }
    }

    public class GetICRandXY : GetICInstance
    {
        private HomogeneousRand rand { get; set; }

        public GetICRandXY() { rand = new HomogeneousRand(); }

        public override Double getX(Double xi)
        {
            return xi * (rand.GetRandomFull() * 2.0 - 1.0);
        }
        public override Double getY(Double yi)
        {
            return yi * (rand.GetRandomFull() * 2.0 - 1.0);
        }
        public override Double getZ(Double zi) { return zi; }
    }

    public class GetICRandXZ : GetICInstance
    {
        private HomogeneousRand rand { get; set; }

        public GetICRandXZ() { rand = new HomogeneousRand(); }

        public override Double getX(Double xi)
        {
            return xi * (rand.GetRandomFull() * 2.0 - 1.0);
        }
        public override Double getZ(Double zi)
        {
            return zi * (rand.GetRandomFull() * 2.0 - 1.0);
        }
        public override Double getY(Double yi) { return yi; }
    }

    public class GetICRandYZ : GetICInstance
    {
        private HomogeneousRand rand { get; set; }

        public GetICRandYZ() { rand = new HomogeneousRand(); }

        public override Double getY(Double yi)
        {
            return yi * (rand.GetRandomFull() * 2.0 - 1.0);
        }
        public override Double getZ(Double zi)
        {
            return zi * (rand.GetRandomFull() * 2.0 - 1.0);
        }
        public override Double getX(Double xi) { return xi; }
    }

    public class GetICRandXYZ : GetICInstance
    {
        private HomogeneousRand rand { get; set; }

        public GetICRandXYZ() { rand = new HomogeneousRand(); }

        public override Double getX(Double xi)
        {
            return xi * (rand.GetRandomFull() * 2.0 - 1.0);
        }
        public override Double getY(Double yi)
        {
            return yi * (rand.GetRandomFull() * 2.0 - 1.0);
        }
        public override Double getZ(Double zi)
        {
            return zi * (rand.GetRandomFull() * 2.0 - 1.0);
        }
    }
}

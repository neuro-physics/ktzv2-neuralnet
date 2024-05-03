using System;

namespace KTzV2.Maths.Matrices.AdjacencyMatrix
{
    public enum AdjacencyMatrixType
    {
        /// <summary>
        /// generates a regular lattice with linear geometry, 2 neighbours and periodic boundary condition
        /// </summary>
        LinearLatticePeriodicBC,

        /// <summary>
        /// generates a regular lattice with linear geometry, 2 neighbours and free boundary condition
        /// </summary>
        LinearLatticeFreeBC,

        /// <summary>
        /// generates a regular lattice with square geometry, 4 neighbours and periodic boundary condition
        /// </summary>
        SquareLatticePeriodicBC,

        /// <summary>
        /// generates a regular lattice with square geometry, 4 neighbours and free boundary condition
        /// </summary>
        SquareLatticeFreeBC,

        /// <summary>
        /// generates a mean-field graph (every element is connected to all the others -- networkX)
        /// </summary>
        CompleteGraph,

        /// <summary>
        /// generates a random graph (networkX)
        /// </summary>
        RandomGraph,

        /// <summary>
        /// generates a Barabasi-Albert Scale-free graph (networkX)
        /// </summary>
        BarabasiAlbertGraph,

        /// <summary>
        /// generates a Watts-Strogatz Small-World graph (it's not guaranteed that all elements are connected in one tree... there may be lonely clusters - poor of them -- networkX)
        /// </summary>
        WattsStrogatzGraph,

        /// <summary>
        /// generates a Watts-Strogatz Small-World graph with every element connected in a single tree... The process can throw an exception if the amount of tries is exceeded (which is 100 -- networkX)
        /// </summary>
        ConnectedWattsStrogatzGraph,

        /// <summary>
        /// generates a cubic lattic with periodic boundary condition - networkX
        /// </summary>
        CubicLatticePeriodicBC,

        /// <summary>
        /// generates a cubic lattic with free boundary condition - networkX
        /// </summary>
        CubicLatticeFreeBC,

        /// <summary>
        /// specifies that the adjacency matrix should be read from file
        /// </summary>
        FromFile
    }

    public class FGetAdjacencyMatrix
    {
        /// <summary>
        /// Get the specified adjacency matrix...
        /// </summary>
        /// <param name="type">the type of the desired matrix</param>
        /// <param name="nElems">total number of elements in the network</param>
        /// <param name="numOfNeighbours">number of neighbours of each element (for a Watts-Strogatz, it's the initial number of neighbours of each element)</param>
        /// <param name="numOfEdgesForNewElem">the "m" parameter for a Barabasi-Albert Graph (num of edges for each new attached node) - only important for Barabasi-Albert Graph</param>
        /// <param name="rewiringProb">the probability to replace a connection with another one (only important to create Watts-Strogatz graph)</param>
        /// <param name="netFileName">name of the file with the adjacency matrix, in case of FromFile selected</param>
        /// <returns>the desired adjacency matrix constructor</returns>
        public static IAdjacencyMatrix<Double> GetMatrixFor(AdjacencyMatrixType type, Int32 nElems, Int32[] L, Int32 numOfNeighbours, Int32 numOfEdgesForNewElem, Double rewiringProb, Boolean isDirected, String netFileName)
        {
            if (type == AdjacencyMatrixType.LinearLatticeFreeBC)
            {
                return new AdjacencyMatrix1DFree(nElems, numOfNeighbours, isDirected);
            }
            else if (type == AdjacencyMatrixType.LinearLatticePeriodicBC)
            {
                return new AdjacencyMatrix1DPeriodic(nElems, numOfNeighbours, isDirected);
            }
            else if (type == AdjacencyMatrixType.SquareLatticeFreeBC)
            {
                if (L[0] == L[1])
                {
                    return new AdjacencyMatrix2DFree(nElems, numOfNeighbours, isDirected);
                }
                else
                {
                    return new AdjacencyMatrixGridGraph(new Int32[] { L[0], L[1] }, false, isDirected);
                }
            }
            else if (type == AdjacencyMatrixType.SquareLatticePeriodicBC)
            {
                if (L[0] == L[1])
                {
                    return new AdjacencyMatrix2DPeriodic(nElems, numOfNeighbours, isDirected);
                }
                else
                {
                    return new AdjacencyMatrixGridGraph(new Int32[] { L[0], L[1] }, true, isDirected);
                }
            }
            else if (type == AdjacencyMatrixType.CompleteGraph)
            {
                return new AdjacencyMatrixCompleteGraph(nElems, isDirected);
            }
            else if (type == AdjacencyMatrixType.RandomGraph)
            {
                return new AdjacencyMatrixRandomGraph(nElems, rewiringProb, isDirected);
            }
            else if (type == AdjacencyMatrixType.BarabasiAlbertGraph)
            {
                return new AdjacencyMatrixBarabasiAlbert(nElems, numOfEdgesForNewElem, isDirected);
            }
            else if (type == AdjacencyMatrixType.WattsStrogatzGraph)
            {
                return new AdjacencyMatrixWattsStrogatz(nElems, numOfNeighbours, rewiringProb, isDirected);
            }
            else if (type == AdjacencyMatrixType.ConnectedWattsStrogatzGraph)
            {
                return new AdjacencyMatrixWattsStrogatzConn(nElems, numOfNeighbours, rewiringProb, isDirected);
            }
            else if (type == AdjacencyMatrixType.CubicLatticeFreeBC)
            {
                return new AdjacencyMatrixGridGraph(L, false, isDirected);
            }
            else if (type == AdjacencyMatrixType.CubicLatticePeriodicBC)
            {
                return new AdjacencyMatrixGridGraph(L, true, isDirected);
            }
            else if (type == AdjacencyMatrixType.FromFile)
            {
                if (netFileName == "") throw new ArgumentNullException("No file specified for adjacency matrix");
                return new AdjacencyMatrixFromFile(netFileName, isDirected);
            }
            throw new ArgumentException("The desired matrix is not available!");
        }
    }

    public interface IAdjacencyMatrix<T> where T : struct
    {
        /// <summary>
        /// number of elements
        /// </summary>
        Int32 NElems { get; }

        /// <summary>
        /// builds and gets the matrix (sparse matrix)
        /// </summary>
        /// <returns>a sparse matrix of type T containing the adjacency matrix</returns>
        SparseMatrix<T> BuildAndGetMatrix();

        /// <summary>
        /// utilized to create directed graphs
        /// </summary>
        SparseMatrix<T> DeleteLowerDiagonalElem(SparseMatrix<T> matrix);
    }

    public sealed class AdjacencyMatrixGridGraph : AdjacencyMatrixFromNetworkX
    {
        /// <summary>
        /// Array containing the amount of elements on each dimension; dim of the graph == dim of this array
        /// </summary>
        public Int32[] L { get; set; }
        /// <summary>
        /// is it a periodic graph
        /// </summary>
        public Boolean isPeriodic { get; set; }
        /// <summary>
        /// constructor
        /// </summary>
        /// <param name="L">Array containing the amount of elements on each dimension; dim of the graph == dim of this array</param>
        /// <param name="isPeriodic">is the graph periodic?</param>
        public AdjacencyMatrixGridGraph(Int32[] L, Boolean isPeriodic, Boolean isDirected)
            : base("grid_graph", isDirected)
        {
            this.L = L;
            this.isPeriodic = isPeriodic;
            this.RunBuildNetProcess();
            base.Initialize(this.FileName);
        }

        /// <summary>
        /// generates the method call according to the chosen GraphMethod and to the NetworkX Docs
        /// </summary>
        /// <returns>a string containing the NetworkX method call</returns>
        protected override String GetNetworkXGeneratorMethodCall()
        {
            String s = "";
            Int32 i = 0;
            while (i < L.Length)
            {
                s += L[i].ToString() + ",";
                i++;
            }
            s = s.Remove(s.Length - 1);
            return String.Format("G=nx.{0}(dim=[{1}],periodic={2})", this.GraphMethod, s, this.IsPeriodic());
        }

        /// <summary>
        /// the python command line to get the adjacency matrix of the selected method from networkx
        /// </summary>
        /// <returns>the python command line to get the adjacency matrix of the selected method from networkx</returns>
        protected override string GetNetworkXAdjMatrix()
        {
            return "nx.adjacency_matrix(G,nodelist=sorted(G.nodes())).todense()";
        }

        /// <summary>
        /// returns a python boolean type in string format to ajust this graph as periodic or not
        /// </summary>
        /// <returns>a string containing a boolean type of the python language conatining the periodicity characteristic of this graph</returns>
        private String IsPeriodic()
        {
            String s = this.isPeriodic.ToString();
            return s[0].ToString().ToUpper() + s.Substring(1);
        }
    }

    public sealed class AdjacencyMatrixCompleteGraph : AdjacencyMatrixFromNetworkX
    {
        /// <summary>
        /// constructor - creates the temporary file with the graph's adj matrix
        /// </summary>
        /// <param name="numOfElems">number of elements in the graph</param>
        public AdjacencyMatrixCompleteGraph(Int32 numOfElems, Boolean isDirected)
            : base("complete_graph", isDirected)
        {
            // setting parameters
            this.NElems = numOfElems;
            this.RunBuildNetProcess();
            base.Initialize(this.FileName);
        }
    }

    public sealed class AdjacencyMatrixWattsStrogatzConn : AdjacencyMatrixFromNetworkX
    {
        /// <summary>
        /// the initial number of neighbours of each element
        /// </summary>
        private Int32 k { get; set; }

        /// <summary>
        /// the rewiring probability
        /// </summary>
        private Double p { get; set; }

        public AdjacencyMatrixWattsStrogatzConn(Int32 numOfElems, Int32 numOfNeighbours, Double rewiringProb, Boolean isDirected)
            : base("connected_watts_strogatz_graph", isDirected)
        {
            // setting parameters
            this.NElems = numOfElems;
            this.k = numOfNeighbours;
            this.p = rewiringProb;
            this.RunBuildNetProcess();
            base.Initialize(this.FileName);
        }

        /// <summary>
        /// generates the method call according to the chosen GraphMethod and to the NetworkX Docs
        /// </summary>
        /// <returns>a string containing the NetworkX method call</returns>
        protected override String GetNetworkXGeneratorMethodCall()
        {
            return "G=nx." + this.GraphMethod + "(" + this.NElems + "," + this.k + "," + this.p + ")";
        }
    }

    public sealed class AdjacencyMatrixWattsStrogatz : AdjacencyMatrixFromNetworkX
    {
        /// <summary>
        /// the initial number of neighbours of each element
        /// </summary>
        private Int32 k { get; set; }

        /// <summary>
        /// the rewiring probability
        /// </summary>
        private Double p { get; set; }

        public AdjacencyMatrixWattsStrogatz(Int32 numOfElems, Int32 numOfNeighbours, Double rewiringProb, Boolean isDirected)
            : base("watts_strogatz_graph", isDirected)
        {
            // setting parameters
            this.NElems = numOfElems;
            this.k = numOfNeighbours;
            this.p = rewiringProb;
            this.RunBuildNetProcess();
            base.Initialize(this.FileName);
        }

        /// <summary>
        /// generates the method call according to the chosen GraphMethod and to the NetworkX Docs
        /// </summary>
        /// <returns>a string containing the NetworkX method call</returns>
        protected override String GetNetworkXGeneratorMethodCall()
        {
            return "G=nx." + this.GraphMethod + "(" + this.NElems + "," + this.k + "," + this.p + ")";
        }
    }

    public sealed class AdjacencyMatrixRandomGraph : AdjacencyMatrixFromNetworkX
    {
        /// <summary>
        /// the initial number of neighbours of each element
        /// </summary>
        private Int32 k { get; set; }

        /// <summary>
        /// the rewiring probability
        /// </summary>
        private Double p { get; set; }

        public AdjacencyMatrixRandomGraph(Int32 numOfElems, Double pOfEdgeCreation, Boolean isDirected)
            : base("fast_gnp_random_graph", isDirected)
        {
            // setting parameters
            this.NElems = numOfElems;
            this.p = pOfEdgeCreation;
            this.RunBuildNetProcess();
            base.Initialize(this.FileName);
        }

        /// <summary>
        /// generates the method call according to the chosen GraphMethod and to the NetworkX Docs
        /// </summary>
        /// <returns>a string containing the NetworkX method call</returns>
        protected override String GetNetworkXGeneratorMethodCall()
        {
            return "G=nx." + this.GraphMethod + "(" + this.NElems + "," + this.p + ", directed=" + (this.IsDirected? "True" : "False") + ")";
        }
    }

    public sealed class AdjacencyMatrixBarabasiAlbert : AdjacencyMatrixFromNetworkX
    {
        /// <summary>
        /// number of edges that a new attached element will have initially
        /// </summary>
        private Int32 m { get; set; }

        /// <summary>
        /// constructor - creates the temporary file with the graph's adj matrix
        /// </summary>
        /// <param name="numOfElems">number of elements in the graph</param>
        /// <param name="numOfEdgesForNewElem">initial number of connections that a newly attached element will have</param>
        public AdjacencyMatrixBarabasiAlbert(Int32 numOfElems, Int32 numOfEdgesForNewElem, Boolean isDirected)
            : base("barabasi_albert_graph", isDirected)
        {
            // setting parameters
            this.NElems = numOfElems;
            this.m = numOfEdgesForNewElem;
            this.RunBuildNetProcess();
            base.Initialize(this.FileName);
        }

        /// <summary>
        /// generates the method call according to the chosen GraphMethod and to the NetworkX Docs
        /// </summary>
        /// <returns>a string containing the NetworkX method call</returns>
        protected override String GetNetworkXGeneratorMethodCall()
        {
            return "G=nx." + this.GraphMethod + "(" + this.NElems + "," + this.m + ")";
        }
    }

    public abstract class AdjacencyMatrixFromNetworkX : AdjacencyMatrixFromFile
    {
        /// <summary>
        /// networkX method used to build the graph
        /// </summary>
        protected String GraphMethod { get; set; }

        /// <summary>
        /// process start info for running python
        /// </summary>
        protected System.Diagnostics.ProcessStartInfo ProcessInfo { get; set; }

        /// <summary>
        /// constructor of the base class
        /// </summary>
        /// <param name="graphMethod">the graph generator method as specified by the NetworkX Documentation</param>
        public AdjacencyMatrixFromNetworkX(String graphMethod, Boolean isDirected)
            : base(isDirected)
        {
            this.GraphMethod = graphMethod;
            this.FileName = this.GraphMethod + ".tmp";
            this.FileName = KTzV2.Data.Header.KTzHeader.CheckAndGetFileName(this.FileName);
        }

        /// <summary>
        /// generates a process start info for this graph
        /// </summary>
        /// <returns></returns>
        protected System.Diagnostics.ProcessStartInfo GetProcessStartInfo(String pythonArgs)
        {
            // creating a process startinfo
            this.ProcessInfo = new System.Diagnostics.ProcessStartInfo();
            this.ProcessInfo.RedirectStandardOutput = true;
            this.ProcessInfo.RedirectStandardError = true;
            this.ProcessInfo.UseShellExecute = false;
            this.ProcessInfo.CreateNoWindow = true;
            this.ProcessInfo.WindowStyle = System.Diagnostics.ProcessWindowStyle.Hidden;
            if (Environment.OSVersion.ToString().Contains("Windows"))
            {
                // if we are in a windows environment, we should run another cmd
                this.ProcessInfo.FileName = "cmd";
                this.ProcessInfo.Arguments = "/c python " + pythonArgs;
            }
            else
            {
                // otherwise, we should just run python
                this.ProcessInfo.FileName = "python";
                this.ProcessInfo.Arguments = pythonArgs;
            }
            return this.ProcessInfo;
        }

        /// <summary>
        /// generates the method call according to the chosen GraphMethod and to the NetworkX Docs
        /// </summary>
        /// <returns>a string containing the NetworkX method call</returns>
        protected virtual String GetNetworkXGeneratorMethodCall()
        {
            return String.Format("G=nx.{0}({1})", this.GraphMethod, this.NElems);
        }

        /// <summary>
        /// the python command line to get the adjacency matrix of the selected method from networkx
        /// </summary>
        /// <returns>the python command line to get the adjacency matrix of the selected method from networkx</returns>
        protected virtual String GetNetworkXAdjMatrix()
        {
            return "nx.adjacency_matrix(G).todense()";
        }

        /// <summary>
        /// runs the process with the networkx
        /// </summary>
        protected void RunBuildNetProcess()
        {
            // creating process
            System.Diagnostics.Process proc = new System.Diagnostics.Process();
            String pythonArgs = String.Format("-c \"import networkx as nx; import numpy as np; {0}; np.savetxt('{1}',{2},fmt='%1d')\"", this.GetNetworkXGeneratorMethodCall(), this.FileName, this.GetNetworkXAdjMatrix());
            proc.StartInfo = this.GetProcessStartInfo(pythonArgs);
            try
            {
                proc.Start();
                proc.WaitForExit();
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
            }
            String result = proc.StandardError.ReadToEnd();
            if (result.Contains("Traceback"))
            {
                Console.WriteLine(result);
                throw new Exception("An error occurred during the graph generation process...");
            }
            //if (System.IO.File.Exists(this.FileName))
            //    System.IO.File.Delete(this.FileName);
        }

        public override SparseMatrix<Double> BuildAndGetMatrix()
        {
            if (!System.IO.File.Exists(this.FileName))
                this.RunBuildNetProcess();
            SparseMatrix<Double> adj = base.BuildAndGetMatrix();
            if (System.IO.File.Exists(this.FileName))
                System.IO.File.Delete(this.FileName);
            return adj;
        }
    }

    public class AdjacencyMatrixFromFile : IAdjacencyMatrix<Double>
    {
        /// <summary>
        /// number of elements
        /// </summary>
        public Int32 NElems { get; protected set; }

        protected Boolean IsDirected { get; set; }

        /// <summary>
        /// the temporary filename where will be stored the adjacency matrix
        /// </summary>
        protected String FileName { get; set; }

        /// <summary>
        /// constructor of the base class
        /// </summary>
        /// <param name="fileName">the name of the file with the adjacency matrix</param>
        public AdjacencyMatrixFromFile(String fileName, Boolean isDirected)
        {
            this.IsDirected = isDirected;
            this.Initialize(fileName);
        }

        public AdjacencyMatrixFromFile(Boolean isDirected)
        {
            this.IsDirected = isDirected;
        }

        /// <summary>
        /// empty constructor
        /// </summary>
        protected AdjacencyMatrixFromFile()
        {
            this.NElems = -1;
        }

        /// <summary>
        /// initializes this adjacency matrix if it has been inherited
        /// </summary>
        /// <param name="fileName">the name of the file with the adjacency matrix</param>
        protected void Initialize(String fileName)
        {
            this.FileName = fileName;
            try
            {
                this.SetNElems();
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
            }
        }

        /// <summary>
        /// adjusts the amount of elements within this adjacency matrix by counting columns of the first line of the specified file
        /// </summary>
        private void SetNElems()
        {
            System.IO.FileStream fs;
            System.IO.StreamReader sr;
            fs = new System.IO.FileStream(this.FileName, System.IO.FileMode.Open, System.IO.FileAccess.Read);
            sr = new System.IO.StreamReader(fs);
            this.NElems = sr.ReadLine().Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries).Length;
            sr.Close();
            fs.Close();
        }

        /// <summary>
        /// builds and gets the matrix (sparse matrix)
        /// </summary>
        /// <returns>a sparse matrix of type T containing the adjacency matrix</returns>
        public virtual SparseMatrix<Double> BuildAndGetMatrix()
        {
            if (this.NElems == -1)
                throw new ArgumentException("AdjacencyMatrixFromFile has not been initialized properly!");

            SparseMatrix<Double> adj = new SparseMatrix<Double>(this.NElems, this.NElems);

            try
            {
                // openning temporary file for reading
                System.IO.FileStream fs;
                System.IO.StreamReader sr;
                fs = new System.IO.FileStream(this.FileName, System.IO.FileMode.Open, System.IO.FileAccess.Read);
                sr = new System.IO.StreamReader(fs);

                String[] col;
                Int32 j, n;
                Double val;
                Int32 i = 0;
                while (!sr.EndOfStream)
                {
                    col = sr.ReadLine().Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                    n = col.Length;
                    j = 0;
                    while (j < n)
                    {
                        val = Convert.ToDouble(col[j]);
                        if (val != 0.0D)
                        {
                            adj.AddElementWithoutCheckingExistence(i, j, val);
                        }
                        j++;
                    }
                    i++;
                }

                sr.Close();
                fs.Close();
            }
            catch (System.IO.IOException e)
            {
                Console.WriteLine(e.Message);
            }
            if (this.IsDirected) // transforms the matrix into a diagonal matrix by simply deleting the lower diagonal
            {
                adj = this.DeleteLowerDiagonalElem(adj);
            }
            return adj;
        }

        public virtual SparseMatrix<Double> DeleteLowerDiagonalElem(SparseMatrix<Double> matrix)
        {
            Int32[] n;
            Double[] elv;
            Int32 i, j;
            i = 0;
            while (i < NElems)
            {
                n = matrix.GetNonEmptyColsInd(i, out elv);
                j = 0;
                while (j < n.Length)
                {
                    if (i > n[j])
                        matrix.RemElement(i, n[j]);
                    j++;
                }
                i++;
            }
            return matrix;
        }
    }

    public abstract class IRegularAdjacencyMatrix : IAdjacencyMatrix<Double>
    {
        //public SparseMatrix<Boolean> matrix { get; protected set; }
        public Int32 NElems { get; protected set; }
        public Int32 NElemsOnARow { get; protected set; }
        public Int32 NNeighbours { get; protected set; }
        public Boolean IsDirected { get; private set; }
        public IRegularAdjacencyMatrix(Int32 nElems, Int32 nElemsOnARow, Int32 nNeighbours, Boolean isDirected)
        {
            this.NNeighbours = nNeighbours;
            this.NElemsOnARow = nElemsOnARow;
            this.NElems = nElems;
            this.IsDirected = isDirected;
        }
        public abstract SparseMatrix<Double> BuildAndGetMatrix();
        public virtual SparseMatrix<Double> DeleteLowerDiagonalElem(SparseMatrix<Double> matrix)
        {
            Int32[] n;
            Double[] elv;
            Int32 i, j;
            i = 0;
            while (i < NElems)
            {
                n = matrix.GetNonEmptyColsInd(i, out elv);
                j = 0;
                while (j < n.Length)
                {
                    if (i > n[j])
                        matrix.RemElement(i, n[j]);
                    j++;
                }
                i++;
            }
            return matrix;
        }
    }

    public sealed class AdjacencyMatrix1DPeriodic : IRegularAdjacencyMatrix
    {
        public AdjacencyMatrix1DPeriodic(Int32 nElems, Int32 nNeighbours, Boolean isDirected)
            : base(nElems, nElems, nNeighbours, isDirected) { }
        public override SparseMatrix<Double> BuildAndGetMatrix()
        {
            SparseMatrix<Double> matrix = new SparseMatrix<Double>(NElems, NElems);
            Int32 i;
            i = 0;
            Int32[] n = new Int32[NNeighbours]; // the index of the neurons connected to each k neuron (i is the row # and j is the column # of the network matrix, which is LxL, here L = nNeuronsOnARow)
            // the first row doesn't connect to n0, the last row doesn't connect to n1, the first column doesn't connect to n2, the last column doesn't connect to n3
            // n0 = (i-1) * L + j
            // n1 = (i+1) * L + j
            // n2 = i * L + (j-1)
            // n3 = i * L + (j+1)
            while (i < NElems)
            {
                if (i != 0)
                {
                    n[0] = i - 1;
                }
                else
                {
                    n[0] = NElems - 1;
                }
                if (i != NElems - 1)
                {
                    n[1] = i + 1;
                }
                else
                {
                    n[1] = 0;
                }
                matrix.AddElementWithoutCheckingExistence(i, n[0], 1.0D);
                matrix.AddElementWithoutCheckingExistence(i, n[1], 1.0D);
                i++;
            }
            if (this.IsDirected) // transforms the matrix into a diagonal matrix by simply deleting the lower diagonal
            {
                matrix = this.DeleteLowerDiagonalElem(matrix);
            }
            return matrix;
        }
    }

    public sealed class AdjacencyMatrix1DFree : IRegularAdjacencyMatrix
    {
        public AdjacencyMatrix1DFree(Int32 nElems, Int32 nNeighbours, Boolean isDirected)
            : base(nElems, nElems, nNeighbours, isDirected) { }
        public override SparseMatrix<Double> BuildAndGetMatrix()
        {
            SparseMatrix<Double> matrix = new SparseMatrix<Double>(NElems, NElems);
            Int32 i;
            i = 0;
            Int32[] n = new Int32[NNeighbours]; // the index of the neurons connected to each k neuron (i is the row # and j is the column # of the network matrix, which is LxL, here L = nNeuronsOnARow)
            // the first row doesn't connect to n0, the last row doesn't connect to n1, the first column doesn't connect to n2, the last column doesn't connect to n3
            // n0 = (i-1) * L + j
            // n1 = (i+1) * L + j
            // n2 = i * L + (j-1)
            // n3 = i * L + (j+1)
            while (i < NElems)
            {
                if (i != 0)
                {
                    n[0] = i - 1; // (i - 1) * nNeurons + j;
                    matrix.AddElementWithoutCheckingExistence(i, n[0], 1.0D);
                }
                if (i != NElems - 1)
                {
                    n[1] = i + 1; // (i + 1) * nNeurons + j;
                    matrix.AddElementWithoutCheckingExistence(i, n[1], 1.0D);
                }
                i++;
            }
            if (this.IsDirected) // transforms the matrix into a diagonal matrix by simply deleting the lower diagonal
            {
                matrix = this.DeleteLowerDiagonalElem(matrix);
            }
            return matrix;
        }
    }

    public sealed class AdjacencyMatrix2DPeriodic : IRegularAdjacencyMatrix
    {
        public AdjacencyMatrix2DPeriodic(Int32 nElems, Int32 nNeighbours, Boolean isDirected)
            : base(nElems, nElems, nNeighbours, isDirected)
        {
            Double nElemsSqrt = Math.Sqrt(nElems);
            if (Math.Floor(nElemsSqrt) == nElemsSqrt)
            {
                this.NElemsOnARow = (Int32)nElemsSqrt;
            }
            else
            {
                throw new ArgumentOutOfRangeException("The specified number of elements does not have an exact sqrt, so it cannot form a sqr net");
            }
        }

        public override SparseMatrix<Double> BuildAndGetMatrix()
        {
            SparseMatrix<Double> matrix = new SparseMatrix<Double>(NElems, NElems);
            Int32 i, j, k, m;
            i = 0;
            Int32[] n = new Int32[NNeighbours]; // the index of the neurons connected to each k neuron (i is the row # and j is the column # of the network matrix, which is LxL, here L = nNeuronsOnARow)
            // the first row doesn't connect to n0, the last row doesn't connect to n1, the first column doesn't connect to n2, the last column doesn't connect to n3
            // n0 = (i-1) * L + j
            // n1 = (i+1) * L + j
            // n2 = i * L + (j-1)
            // n3 = i * L + (j+1)
            while (i < NElemsOnARow)
            {
                j = 0;
                while (j < NElemsOnARow)
                {
                    // the index of the neuron at site i,j on the network
                    k = i * NElemsOnARow + j;

                    // the indeces of the neurons around it
                    n[0] = (i - 1) * NElemsOnARow + j;
                    n[1] = (i + 1) * NElemsOnARow + j;
                    n[2] = i * NElemsOnARow + (j - 1);
                    n[3] = i * NElemsOnARow + (j + 1);

                    // the matrix should be symmetric (Aij = Aji, because every connection is a two-way connection)
                    m = 0;
                    while (m < NNeighbours)
                    {
                        if ((i == 0) && (m == 0)) // first row connects to last row
                        {
                            n[0] = (NElemsOnARow - 1) * NElemsOnARow + j; // i = nNeuronsOnARow - 1
                        }
                        if ((i == (NElemsOnARow - 1)) && (m == 1)) // last row connects to the first row
                        {
                            n[1] = j; // i = 0
                        }
                        if ((j == 0) && (m == 2)) // first column connects to the last column
                        {
                            n[2] = i * NElemsOnARow + (NElemsOnARow - 1); // j = nNeuronsOnARow - 1
                        }
                        if ((j == (NElemsOnARow - 1)) && (m == 3)) // last column connects to the first column
                        {
                            n[3] = i * NElemsOnARow; // j = 0
                        }
                        matrix.AddElementWithoutCheckingExistence(k, n[m], 1.0D);
                        //AMatrix.addElement(n[m], k, true);
                        m++;
                    }
                    j++;
                }
                i++;
            }
            if (this.IsDirected) // transforms the matrix into a diagonal matrix by simply deleting the lower diagonal
            {
                matrix = this.DeleteLowerDiagonalElem(matrix);
            }
            return matrix;
        }
    }

    public sealed class AdjacencyMatrix2DFree : IRegularAdjacencyMatrix
    {
        public AdjacencyMatrix2DFree(Int32 nElems, Int32 nNeighbours, Boolean isDirected)
            : base(nElems, nElems, nNeighbours, isDirected)
        {
            Double nElemsSqrt = Math.Sqrt(nElems);
            if (Math.Floor(nElemsSqrt) == nElemsSqrt)
            {
                this.NElemsOnARow = (Int32)nElemsSqrt;
            }
            else
            {
                throw new ArgumentOutOfRangeException("The specified number of elements does not have an exact sqrt, so it cannot form a sqr net");
            }
        }
        public override SparseMatrix<Double> BuildAndGetMatrix()
        {
            SparseMatrix<Double> matrix = new SparseMatrix<Double>(NElems, NElems);
            Int32 i, j, k, m;
            i = 0;
            Int32[] n = new Int32[NNeighbours]; // the index of the neurons connected to each k neuron (i is the row # and j is the column # of the network matrix, which is LxL, here L = nNeuronsOnARow)
            // the first row doesn't connect to n0, the last row doesn't connect to n1, the first column doesn't connect to n2, the last column doesn't connect to n3
            // n0 = (i-1) * L + j
            // n1 = (i+1) * L + j
            // n2 = i * L + (j-1)
            // n3 = i * L + (j+1)
            while (i < NElemsOnARow)
            {
                j = 0;
                while (j < NElemsOnARow)
                {
                    // the index of the neuron at site i,j on the network
                    k = i * NElemsOnARow + j;

                    // the indeces of the neurons around it
                    n[0] = (i - 1) * NElemsOnARow + j;
                    n[1] = (i + 1) * NElemsOnARow + j;
                    n[2] = i * NElemsOnARow + (j - 1);
                    n[3] = i * NElemsOnARow + (j + 1);

                    // the matrix should be symmetric (Aij = Aji, because every connection is a two-way connection)
                    m = 0;
                    while (m < NNeighbours)
                    {
                        if ((i == 0) && (m == 0)) // first row doesn't connect to n0
                        {
                            m++;
                            continue;
                        }
                        if ((i == (NElemsOnARow - 1)) && (m == 1)) // last row doesn't connect to n1
                        {
                            m++;
                            continue;
                        }
                        if ((j == 0) && (m == 2)) // first column doesn't connect to n2
                        {
                            m++;
                            continue;
                        }
                        if ((j == (NElemsOnARow - 1)) && (m == 3)) // last column doesn't connect to n3
                        {
                            m++;
                            continue;
                        }
                        matrix.AddElementWithoutCheckingExistence(k, n[m], 1.0D);
                        //AMatrix.addElement(n[m], k, true);
                        m++;
                    }
                    j++;
                }
                i++;
            }
            if (this.IsDirected) // transforms the matrix into a diagonal matrix by simply deleting the lower diagonal
            {
                matrix = this.DeleteLowerDiagonalElem(matrix);
            }
            return matrix;
        }
    }
}

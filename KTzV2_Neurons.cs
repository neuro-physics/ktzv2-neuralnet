using System;
using System.Collections.Generic;
using KTzV2.Synapses;

namespace KTzV2.Neurons
{
    public enum NeuronType
    {
        KTH,
        KTHLog,
        KTz,
        KTzLog,
        KTz2Tanh,
        KTzMF,
        KTzLogMF,
        SIElement,
        GLNeuron,
        GLNeuronLHG,
        GLNeuronLHG1Par,
        GLNeuronFacilitation,
        GLNeuronGammaStochastic
    }

    public class NeuronParam
    {
        public Int32 ind { get; private set; }
        public Double K { get; private set; }
        public Double T { get; private set; }
        public Double d { get; private set; }
        public Double l { get; private set; }
        public Double xR { get; private set; }
        public Double H { get; private set; }
        public Double Q { get; private set; }
        public Double x0 { get; private set; }
        public Double y0 { get; private set; }
        public Double z0 { get; private set; }
        public Double Theta { get; private set; }
        public Double VB { get; private set; }
        public Double VR { get; private set; }
        public Double VT { get; private set; }
        public Double mu { get; private set; }
        public Double tau1 { get; private set; }
        public Double Gamma1 { get; private set; }
        public Double tau0 { get; private set; }
        public Double Gamma0 { get; private set; }
        public Double Gamma_init { get; private set; }
        public Double V0 { get; private set; }
        public Double pS { get; private set; }
        public Double pR { get; private set; }
        public Boolean normalizeInput { get; private set; }
        public UInt64? seed { get; private set; }

        /// <param name="ind">index of the neuron inside a network</param>
        /// <param name="K">parameter K</param>
        /// <param name="T">parameter T</param>
        /// <param name="d">parameter delta</param>
        /// <param name="l">parameter lambda</param>
        /// <param name="xR">parameter xR</param>
        /// <param name="H">parameter H</param>
        /// <param name="Theta">Threshold for SIElement</param>
        /// <param name="x_init">initial value of x</param>
        /// <param name="y_init">initial value of y</param>
        /// <param name="z_init">initial value of z</param>
        public NeuronParam(Int32 ind, Double K, Double T, Double d, Double l, Double xR, Double H, Double x_init, Double y_init, Double z_init, Double Theta,
            Double tau0 = 0, Double Gamma0 = 0, Double tau1 = 0, Double Gamma1 = 0, Double pS = 0, Double pR = 0, Double VR = 0, Double VT = 0, Double VB = 0,
            Double mu = 0, Double V_init = 0, Double Gamma_init = 0, Boolean normalizeInput = false, UInt64? seed = null,Double Q = 0.0D)
        {
            this.ind = ind;
            this.K = K;
            this.T = T;
            this.d = d;
            this.l = l;
            this.xR = xR;
            this.H = H;
            this.Q = Q;
            this.x0 = x_init;
            this.y0 = y_init;
            this.z0 = z_init;
            this.Theta = Theta;
            this.tau0 = tau0;
            this.Gamma0 = Gamma0;
            this.tau1 = tau1;
            this.Gamma1 = Gamma1;
            this.pS = pS;
            this.pR = pR;
            this.VR = VR;
            this.VT = VT;
            this.VB = VB;
            this.mu = mu;
            this.V0 = V_init;
            this.Gamma_init = Gamma_init;
            this.normalizeInput = normalizeInput;
            this.seed = seed;
        }
    }

    public static class NeuronFactory
    {
        public static INeuron GetKTNeuron(NeuronType neuron, NeuronParam par)
        {
            switch (neuron)
            {
                case NeuronType.KTH:
                    return new KTSimple(par.ind, par.K, par.T, par.H, par.x0, par.y0);
                case NeuronType.KTHLog:
                    return new KTLogSimple(par.ind, par.K, par.T, par.H, par.x0, par.y0);
                case NeuronType.KTz:
                    return new KTzNeuron(par.ind, par.K, par.T, par.d, par.l, par.xR, par.x0, par.y0, par.z0);
                case NeuronType.KTz2Tanh:
                    return new KTz2Tanh(par.ind, par.K, par.T, par.d, par.l, par.xR, par.H, par.Q, par.x0, par.y0, par.z0);
                case NeuronType.KTzLog:
                    return new KTzLogNeuron(par.ind, par.K, par.T, par.d, par.l, par.xR, par.x0, par.y0, par.z0);
                case NeuronType.KTzMF:
                    return new KTzNeuronAvgInp(par.ind, par.K, par.T, par.d, par.l, par.xR, par.x0, par.y0, par.z0);
                case NeuronType.KTzLogMF:
                    return new KTzLogNeuronAvgInp(par.ind, par.K, par.T, par.d, par.l, par.xR, par.x0, par.y0, par.z0);
                case NeuronType.SIElement:
                    return new SIElement(par.ind, par.Theta, par.x0);
                case NeuronType.GLNeuron:
                    return new GLNeuronLogistic(par.ind, par.VR, par.VT, par.VB, par.mu, par.x0, par.V0, par.Gamma_init, par.normalizeInput, par.seed);
                case NeuronType.GLNeuronLHG:
                    return new GLNeuronAdaptable(par.ind, 0.0D, 1.0D / par.tau0, 0.0D, par.tau1, par.Gamma1, par.VR, par.VT, par.VB, par.mu, par.x0, par.V0, par.Gamma_init, par.normalizeInput, par.seed);
                case NeuronType.GLNeuronLHG1Par:
                    return new GLNeuronAdaptable(par.ind, 0.0D, 1.0D, 0.0D, -par.tau1, 0.0D, par.VR, par.VT, par.VB, par.mu, par.x0, par.V0, par.Gamma_init, par.normalizeInput, par.seed);
                case NeuronType.GLNeuronFacilitation:
                    return new GLNeuronAdaptable(par.ind, 1.0D, par.tau0, par.Gamma0, par.tau1, par.Gamma1, par.VR, par.VT, par.VB, par.mu, par.x0, par.V0, par.Gamma_init, par.normalizeInput, par.seed);
                case NeuronType.GLNeuronGammaStochastic:
                    return new GLNeuronGammaStochastic(par.ind, par.pR, par.Gamma0, par.pS, par.Gamma1, par.VR, par.VT, par.VB, par.mu, par.x0, par.V0, par.Gamma_init, par.normalizeInput, par.seed);
                default:
                    throw new ArgumentOutOfRangeException("Unrecognized neuron type! " + neuron.ToString());
            }
        }
    }

    public class GLNeuronAdaptable : GLNeuronLogistic
    {
        /// <summary>
        /// growth time constant of Gamma
        /// </summary>
        public Double iTau1 { get; protected set; }
        /// <summary>
        /// Gamma superior asymptotic value
        /// </summary>
        public Double Gamma1 { get; protected set; }
        /// <summary>
        /// decay time constant of Gamma
        /// </summary>
        public Double iTau0 { get; protected set; }
        /// <summary>
        /// Gamma inferior asymptotic value
        /// </summary>
        public Double Gamma0 { get; protected set; }
        /// <summary>
        /// dynamic selector (delta == 0 -> LHG-like, delta == 1 -> Facilitation)
        /// </summary>
        public Double delta { get; protected set; }

        public GLNeuronAdaptable(Int32 ind, Double delta, Double tau0, Double Gamma0, Double tau1, Double Gamma1, Double VR, Double VT, Double VB, Double mu, Double x_init, Double V_init, Double Gamma_init, Boolean normalizeInput = false, UInt64? seed = null)
            : base(ind, VR, VT, VB, mu, x_init, V_init, Gamma_init, normalizeInput, seed)
        {
            this.delta = delta;
            this.iTau1 = 1 / tau1;
            this.iTau0 = 1 / tau0;
            this.Gamma1 = Gamma1;
            this.Gamma0 = Gamma0;
        }
        /// <summary>
        /// evaluates all this neuron's variables for one timestep
        /// </summary>
        public override void Evolve()
        {
            this.Evolve(0.0D);
        }

        /// <summary>
        /// evaluates this neuron's dynamical variables with external stimulus I
        /// </summary>
        /// <param name="I">the external stimulus on the membrane potential</param>
        public override void Evolve(Double I)
        {
            this.Gamma = this.Gamma + this.iTau1 * (this.Gamma1 - this.Gamma) * (this.delta * this.x + 1 - this.delta) - this.iTau0 * (this.Gamma - this.Gamma0) * ((1 - this.delta) * this.x + this.delta);
            base.Evolve(I);
        }
    }

    public class GLNeuronGammaStochastic : GLNeuronLogistic
    {
        /// <summary>
        /// seizure probability (p of setting Gamma == Gamma1)
        /// </summary>
        public Double pS { get; protected set; }
        /// <summary>
        /// Gamma superior asymptotic value
        /// </summary>
        public Double Gamma1 { get; protected set; }
        /// <summary>
        /// recovery probability (p of setting Gamma == Gamma0)
        /// </summary>
        public Double pR { get; protected set; }
        /// <summary>
        /// Gamma inferior asymptotic value
        /// </summary>
        public Double Gamma0 { get; protected set; }
        /// <summary>
        /// random number to decide state of Gamma
        /// </summary>
        private KTzV2.Maths.Random.HomogeneousRand randG { get; set; }

        public GLNeuronGammaStochastic(Int32 ind, Double pR, Double Gamma0, Double pS, Double Gamma1, Double VR, Double VT, Double VB, Double mu, Double x_init, Double V_init, Double Gamma_init, Boolean normalizeInput = false, UInt64? seed = null)
            : base(ind, VR, VT, VB, mu, x_init, V_init, Gamma_init, normalizeInput, seed)
        {
            this.pS = pS;
            this.pR = pR;
            this.Gamma1 = Gamma1;
            this.Gamma0 = Gamma0;
            UInt64 s = getSpecificSeed(pS, pR, Gamma0, Gamma1);
            if (seed.HasValue)
                this.randG = new KTzV2.Maths.Random.HomogeneousRand(seed.Value + s + (UInt64)ind);
            else
                this.randG = new KTzV2.Maths.Random.HomogeneousRand(Convert.ToUInt64(DateTime.Now.Ticks) + s + (UInt64)ind);
        }
        private UInt64 getSpecificSeed(Double a, Double b, Double c, Double d)
        {
            String s = a.ToString("0.0000000000E+000") + b.ToString("0.0000000000E+000") + c.ToString("0.0000000000E+000") + d.ToString("0.0000000000E+000");
            System.Security.Cryptography.MD5 m = System.Security.Cryptography.MD5.Create();
            m.Initialize();
            m.ComputeHash(System.Text.Encoding.UTF8.GetBytes(s));
            return UInt64.Parse(m.Hash.ToString().Substring(0, 16), System.Globalization.NumberStyles.HexNumber) + UInt64.Parse(m.Hash.ToString().Substring(16), System.Globalization.NumberStyles.HexNumber);
        }
        /// <summary>
        /// evaluates all this neuron's variables for one timestep
        /// </summary>
        public override void Evolve()
        {
            this.Evolve(0.0D);
        }

        /// <summary>
        /// evaluates this neuron's dynamical variables with external stimulus I
        /// </summary>
        /// <param name="I">the external stimulus on the membrane potential</param>
        public override void Evolve(Double I)
        {
            Double rS = this.randG.GetRandomFull();
            if (rS < this.pS) // transition to seizure state
                this.Gamma = this.Gamma1;
            else if ((rS >= this.pS) && (rS < (this.pS+this.pR))) // transition to healthy state
                this.Gamma = this.Gamma0;
            base.Evolve(I);
        }
    }

    /// <summary>
    /// GL Neuron unit using the function phi = (Gamma*V)/(1+Gamma*V)
    /// </summary>
    public class GLNeuronLogistic : StochasticGLElement
    {
        public GLNeuronLogistic(Int32 ind, Double VR, Double VT, Double VB, Double mu, Double x_init, Double V_init, Double Gamma_init, Boolean normalizeInput = false, UInt64? seed = null)
            : base(ind,VR,VT,VB,mu,x_init,V_init,Gamma_init,normalizeInput,seed)
        {
        }

        /// <summary>
        /// Logistic firing probability
        /// </summary>
        /// <param name="V">membrane potential</param>
        /// <returns></returns>
        protected override Double FiringProbability(double V)
        {
            return this.Gamma * (this.V - this.VT) * (this.V > this.VT? 1.0D : 0.0D) / (1.0D + this.Gamma * (this.V - this.VT));
        }

        /// <summary>
        /// evaluates all this neuron's variables for one timestep
        /// </summary>
        public override void Evolve()
        {
            this.SumInputSignal(); // sums up to this.Isyn
            if (this.x == 0)
            {
                this.V = this.mu * (this.V - this.VB) + this.VB + this.Isyn;
            }
            else
            {
                this.V = this.VR;
            }
            this.x = (this.rand.GetRandomFull() < this.FiringProbability(this.V) ? 1.0D : 0.0D);
        }

        /// <summary>
        /// evaluates this neuron's dynamical variables with external stimulus I
        /// </summary>
        /// <param name="I">the external stimulus on the membrane potential</param>
        public override void Evolve(Double I)
        {
            this.SumInputSignal(); // sums up to this.Isyn
            if (this.x == 0)
            {
                this.V = this.mu * (this.V - this.VB) + this.VB + this.Isyn + I;
            }
            else
            {
                this.V = this.VR;
            }
            this.x = (this.rand.GetRandomFull() < this.FiringProbability(this.V) ? 1.0D : 0.0D);
        }
    }

    public abstract class StochasticGLElement : INeuron
    {
        /// <summary>
        /// the index of this neuron (a number that identifies it within a network)
        /// </summary>
        public Int32 Index { get; protected set; }

        /// <summary>
        /// the list of this neuron's neighbours (only the ones that influences this one)
        /// </summary>
        public List<ISynapse> Input { get; private set; }

        /// <summary>
        /// the sum over all neighbours output signals
        /// </summary>
        public Double Isyn { get; private set; }

        /// <summary>
        /// the sum over all neighbours connection weight
        /// </summary>
        public Double TotalInputWeight { get; private set; }

        /// <summary>
        /// if normalizeInput == true, then the input signal is divided by TotalInputWeight
        /// </summary>
        public Boolean normalizeInput { get; private set; }

        /// <summary>
        /// state variable (0 or 1)
        /// </summary>
        public Double x { get; protected set; }

        /// <summary>
        /// membrane potential used to calculate firing probability
        /// </summary>
        public Double V { get; protected set; }

        /// <summary>
        /// Gain of the firing probability (parameter or variable)
        /// </summary>
        public Double Gamma { get; protected set; }

        /// <summary>
        /// Time constant for V integration
        /// </summary>
        public Double mu { get; protected set; }

        /// <summary>
        /// Baseline potential (if no activity, V = VB)
        /// </summary>
        public Double VB { get; protected set; }

        /// <summary>
        /// Reset potential (V = VR if neuron spikes)
        /// </summary>
        public Double VR { get; protected set; }

        /// <summary>
        /// Firing threshold
        /// </summary>
        public Double VT { get; protected set; }

        protected KTzV2.Maths.Random.HomogeneousRand rand { get; private set; }

        /// <summary>
        /// creates a stochastic neuron of the Galves-Locherbach type
        /// </summary>
        /// <param name="ind">index of this neuron in the network</param>
        /// <param name="VR">reset membrane potential</param>
        /// <param name="VT">threshold</param>
        /// <param name="VB">baseline potential</param>
        /// <param name="mu">dissipation constant</param>
        /// <param name="x_init">spike variable Init condition</param>
        /// <param name="V_init">membrane potential initial condition</param>
        /// <param name="Gamma_init">initial value for gain constant</param>
        /// <param name="normalizeInput">chooses to normalize input</param>
        /// <param name="seed">seed for random number generator</param>
        protected StochasticGLElement(Int32 ind, Double VR, Double VT, Double VB, Double mu, Double x_init, Double V_init, Double Gamma_init, Boolean normalizeInput = false, UInt64? seed = null)
        {
            this.Index = ind;
            this.x = x_init;
            this.V = V_init;
            this.Gamma = Gamma_init;
            this.VB = VB;
            this.VR = VR;
            this.VT = VT;
            this.mu = mu;
            this.normalizeInput = normalizeInput;
            this.SetTotalInputWeight();
            if (seed.HasValue)
                this.rand = new KTzV2.Maths.Random.HomogeneousRand(seed.Value + (UInt64)this.Index);
            else
                this.rand = new KTzV2.Maths.Random.HomogeneousRand(Convert.ToUInt64(DateTime.Now.Ticks) + (UInt64)this.Index);
            this.Input = new List<ISynapse>();
        }

        /// <summary>
        /// Phi(V) function that gives the firing probability
        /// </summary>
        /// <param name="V">membrane potential</param>
        /// <returns>firing probability</returns>
        protected abstract Double FiringProbability(Double V);

        private void SetTotalInputWeight()
        {
            this.TotalInputWeight = 1.0D;
            if (this.normalizeInput)
            {
                Int32 i = 0;
                while (i < this.Input.Count)
                {
                    this.TotalInputWeight += this.Input[i].GetCoupling();
                    i++;
                }
            }
        }

        /// <summary>
        /// initializes the neighbours with the specified list
        /// </summary>
        /// <param name="s">the list of input synapses</param>
        public void AddInput(List<ISynapse> s)
        {
            this.Input = s;
            this.SetTotalInputWeight();
        }

        /// <summary>
        /// adds a neighbour to this neuron...
        /// </summary>
        /// <param name="s">the input synapse</param>
        public void AddInput(ISynapse s)
        {
            this.Input.Add(s);
            this.SetTotalInputWeight();
        }

        /// <summary>
        /// sums all the neighbours output signals and stores it at this.Isyn
        /// </summary>
        protected void SumInputSignal()
        {
            Int32 i = 0, n = this.Input.Count;
            this.Isyn = 0.0;
            while (i < n)
            {
                this.Isyn += this.Input[i].GetSignal();
                i++;
            }
        }

        /// <summary>
        /// returns the membrane potential of the neuron
        /// </summary>
        /// <returns>this neuron membrane potential</returns>
        public Double GetMembranePotential()
        {
            return this.x;
        }

        /// <summary>
        /// evaluates all this neuron's variables for one timestep
        /// </summary>
        public abstract void Evolve();

        /// <summary>
        /// evaluates this neuron's dynamical variables with external stimulus I
        /// </summary>
        /// <param name="I">the external stimulus on the membrane potential</param>
        public abstract void Evolve(Double I);

        /// <summary>
        /// resets neuron variables to the values contained in ic array
        /// </summary>
        /// <param name="ic">new values for neurons variables</param>
        public virtual void ResetIC(Double[] ic)
        {
            this.x = ic[0];
        }

        public virtual bool SpikeDetector()
        {
            return this.x == 1.0;
        }
    }

    /// <summary>
    /// implements the neuron of the reference A simple model for excitable and bursting elements
    /// </summary>
    public class KTz2Tanh : KTSimple
    {
        /// <summary>
        /// parameter delta
        /// </summary>
        public Double d { get; private set; }

        /// <summary>
        /// parameter lambda
        /// </summary>
        public Double l { get; private set; }

        /// <summary>
        /// parameter x_R
        /// </summary>
        public Double xR { get; private set; }

        /// <summary>
        /// parameter H external field in Y
        /// </summary>
        public new Double H { get; private set; }

        /// <summary>
        /// parameter Q external field in X
        /// </summary>
        public Double Q { get; private set; }

        /// <summary>
        /// dynamical variable z
        /// </summary>
        public Double z { get; private set; }

        /// <summary>
        /// initial value of the variable z
        /// </summary>
        public Double z_init { get; private set; }

        /// <summary>
        /// KTzLogNeuron constructor
        /// </summary>
        /// <param name="K">parameter K</param>
        /// <param name="T">parameter T</param>
        /// <param name="d">parameter delta</param>
        /// <param name="l">parameter lambda</param>
        /// <param name="xR">parameter xR</param>
        /// <param name="H">parameter H external field in Y</param>
        /// <param name="Q">parameter Q external field in X</param>
        /// <param name="x_init">initial value of x</param>
        /// <param name="y_init">initial value of y</param>
        /// <param name="z_init">initial value of z</param>
        public KTz2Tanh(Int32 ind, Double K, Double T, Double d, Double l, Double xR, Double H, Double Q, Double x_init, Double y_init, Double z_init)
            : base(ind,K,T,H,x_init,y_init)
        {
            this.d = d;
            this.l = l;
            this.xR = xR;
            this.H = H;
            this.Q = Q;
            this.z_init = this.z = z_init;
        }

        /// <summary>
        /// evaluates this neuron for one timestep
        /// </summary>
        public override void Evolve()
        {
            this.SumInputSignal();
            this.x_prev = this.x;
            this.x = Math.Tanh((this.x - this.K * this.y + this.z + this.Q + this.Isyn) / this.T);
            this.y = Math.Tanh((this.x_prev + this.H)/this.T);
            this.z = (1.0 - this.d) * this.z - this.l * (this.x_prev - this.xR);
        }

        /// <summary>
        /// evaluates this neuron one timestep with an external stimulus I
        /// </summary>
        /// <param name="I">external stimulus</param>
        public override void Evolve(Double I)
        {
            this.SumInputSignal();
            this.x_prev = this.x;
            this.x = Math.Tanh((this.x - this.K * this.y + this.z + this.Q + this.Isyn + I) / this.T);
            this.y = Math.Tanh((this.x_prev + this.H) / this.T);
            this.z = (1.0 - this.d) * this.z - this.l * (this.x_prev - this.xR);
        }

        /// <summary>
        /// resets neuron variables to the values contained in ic array
        /// </summary>
        /// <param name="ic">new values for neurons variables</param>
        public override void ResetIC(Double[] ic)
        {
            base.ResetIC(ic);
            this.z_init = ic[2];
            this.z = ic[2];
        }
    }

    /// <summary>
    /// implements the neuron of the reference A simple model for excitable and bursting elements
    /// </summary>
    public class KTzLogNeuron : KTNeuron
    {
        /// <summary>
        /// parameter delta
        /// </summary>
        public Double d { get; private set; }

        /// <summary>
        /// parameter lambda
        /// </summary>
        public Double l { get; private set; }

        /// <summary>
        /// parameter x_R
        /// </summary>
        public Double xR { get; private set; }

        /// <summary>
        /// dynamical variable z
        /// </summary>
        public Double z { get; private set; }

        /// <summary>
        /// initial value of the variable z
        /// </summary>
        public Double z_init { get; private set; }

        /// <summary>
        /// KTzLogNeuron constructor
        /// </summary>
        /// <param name="K">parameter K</param>
        /// <param name="T">parameter T</param>
        /// <param name="d">parameter delta</param>
        /// <param name="l">parameter lambda</param>
        /// <param name="xR">parameter xR</param>
        /// <param name="x_init">initial value of x</param>
        /// <param name="y_init">initial value of y</param>
        /// <param name="z_init">initial value of z</param>
        public KTzLogNeuron(Int32 ind, Double K, Double T, Double d, Double l, Double xR, Double x_init, Double y_init, Double z_init)
            : base(ind, K, T, x_init, y_init)
        {
            this.d = d;
            this.l = l;
            this.xR = xR;
            this.z_init = this.z = z_init;
        }

        private Double logFun(Double u)
        {
            return u / (1 + (u > 0.0D ? u : -u));
        }

        /// <summary>
        /// evaluates this neuron for one timestep
        /// </summary>
        public override void Evolve()
        {
            this.SumInputSignal();
            this.x_prev = this.x;
            this.x = logFun((this.x - this.K * this.y + this.z + this.Isyn) / this.T);
            this.y = this.x_prev;
            this.z = (1.0 - this.d) * this.z - this.l * (this.x_prev - this.xR);
        }

        /// <summary>
        /// evaluates this neuron one timestep with an external stimulus I
        /// </summary>
        /// <param name="I">external stimulus</param>
        public override void Evolve(Double I)
        {
            this.SumInputSignal();
            this.x_prev = this.x;
            this.x = logFun((this.x - this.K * this.y + this.z + this.Isyn + I) / this.T);
            this.y = this.x_prev;
            this.z = (1.0 - this.d) * this.z - this.l * (this.x_prev - this.xR);
        }

        /// <summary>
        /// resets neuron variables to the values contained in ic array
        /// </summary>
        /// <param name="ic">new values for neurons variables</param>
        public override void ResetIC(Double[] ic)
        {
            base.ResetIC(ic);
            this.z_init = ic[2];
            this.z = ic[2];
        }
    }

    /// <summary>
    /// implements the neuron of the reference A simple model for excitable and bursting elements
    /// </summary>
    public class KTzNeuron : KTNeuron
    {
        /// <summary>
        /// parameter delta
        /// </summary>
        public Double d { get; private set; }

        /// <summary>
        /// parameter lambda
        /// </summary>
        public Double l { get; private set; }

        /// <summary>
        /// parameter x_R
        /// </summary>
        public Double xR { get; private set; }

        /// <summary>
        /// dynamical variable z
        /// </summary>
        public Double z { get; private set; }

        /// <summary>
        /// initial value of the variable z
        /// </summary>
        public Double z_init { get; private set; }

        /// <summary>
        /// KTzNeuron constructor
        /// </summary>
        /// <param name="K">parameter K</param>
        /// <param name="T">parameter T</param>
        /// <param name="d">parameter delta</param>
        /// <param name="l">parameter lambda</param>
        /// <param name="xR">parameter xR</param>
        /// <param name="x_init">initial value of x</param>
        /// <param name="y_init">initial value of y</param>
        /// <param name="z_init">initial value of z</param>
        public KTzNeuron(Int32 ind, Double K, Double T, Double d, Double l, Double xR, Double x_init, Double y_init, Double z_init)
            : base(ind, K, T, x_init, y_init)
        {
            this.d = d;
            this.l = l;
            this.xR = xR;
            this.z_init = this.z = z_init;
        }

        /// <summary>
        /// evaluates this neuron for one timestep
        /// </summary>
        public override void Evolve()
        {
            this.SumInputSignal();
            this.x_prev = this.x;
            this.x = Math.Tanh((this.x - this.K * this.y + this.z + this.Isyn) / this.T);
            this.y = this.x_prev;
            this.z = (1.0 - this.d) * this.z - this.l * (this.x_prev - this.xR);
        }

        /// <summary>
        /// evaluates this neuron one timestep with an external stimulus I
        /// </summary>
        /// <param name="I">external stimulus</param>
        public override void Evolve(Double I)
        {
            this.SumInputSignal();
            this.x_prev = this.x;
            this.x = Math.Tanh((this.x - this.K * this.y + this.z + this.Isyn + I) / this.T);
            this.y = this.x_prev;
            this.z = (1.0 - this.d) * this.z - this.l * (this.x_prev - this.xR);
        }

        /// <summary>
        /// resets neuron variables to the values contained in ic array
        /// </summary>
        /// <param name="ic">new values for neurons variables</param>
        public override void ResetIC(Double[] ic)
        {
            base.ResetIC(ic);
            this.z_init = ic[2];
            this.z = ic[2];
        }
    }

    /// <summary>
    /// implements the neuron of the reference A simple model for excitable and bursting elements with average input
    /// </summary>
    public class KTzLogNeuronAvgInp : KTNeuron
    {
        /// <summary>
        /// parameter delta
        /// </summary>
        public Double d { get; private set; }

        /// <summary>
        /// parameter lambda
        /// </summary>
        public Double l { get; private set; }

        /// <summary>
        /// parameter x_R
        /// </summary>
        public Double xR { get; private set; }

        /// <summary>
        /// dynamical variable z
        /// </summary>
        public Double z { get; private set; }

        /// <summary>
        /// initial value of the variable z
        /// </summary>
        public Double z_init { get; private set; }

        /// <summary>
        /// KTzLogNeuronAvgInp constructor
        /// </summary>
        /// <param name="K">parameter K</param>
        /// <param name="T">parameter T</param>
        /// <param name="d">parameter delta</param>
        /// <param name="l">parameter lambda</param>
        /// <param name="xR">parameter xR</param>
        /// <param name="x_init">initial value of x</param>
        /// <param name="y_init">initial value of y</param>
        /// <param name="z_init">initial value of z</param>
        public KTzLogNeuronAvgInp(Int32 ind, Double K, Double T, Double d, Double l, Double xR, Double x_init, Double y_init, Double z_init)
            : base(ind, K, T, x_init, y_init)
        {
            this.d = d;
            this.l = l;
            this.xR = xR;
            this.z_init = this.z = z_init;
        }
        private Double logFun(Double u)
        {
            return u / (1 + (u > 0.0D ? u : -u));
        }

        /// <summary>
        /// evaluates this neuron for one timestep
        /// </summary>
        public override void Evolve()
        {
            this.SumInputSignal();
            this.x_prev = this.x;
            this.x = logFun((this.x - this.K * this.y + this.z + this.Isyn / this.Input.Count) / this.T);
            this.y = this.x_prev;
            this.z = (1.0 - this.d) * this.z - this.l * (this.x_prev - this.xR);
        }

        /// <summary>
        /// evaluates this neuron one timestep with an external stimulus I
        /// </summary>
        /// <param name="I">external stimulus</param>
        public override void Evolve(Double I)
        {
            this.SumInputSignal();
            this.x_prev = this.x;
            this.x = logFun((this.x - this.K * this.y + this.z + I + this.Isyn / (this.Input.Count + 1)) / this.T);
            this.y = this.x_prev;
            this.z = (1.0 - this.d) * this.z - this.l * (this.x_prev - this.xR);
        }

        /// <summary>
        /// resets neuron variables to the values contained in ic array
        /// </summary>
        /// <param name="ic">new values for neurons variables</param>
        public override void ResetIC(Double[] ic)
        {
            base.ResetIC(ic);
            this.z_init = ic[2];
            this.z = ic[2];
        }
    }

    /// <summary>
    /// implements the neuron of the reference A simple model for excitable and bursting elements with average input
    /// </summary>
    public class KTzNeuronAvgInp : KTNeuron
    {
        /// <summary>
        /// parameter delta
        /// </summary>
        public Double d { get; private set; }

        /// <summary>
        /// parameter lambda
        /// </summary>
        public Double l { get; private set; }

        /// <summary>
        /// parameter x_R
        /// </summary>
        public Double xR { get; private set; }

        /// <summary>
        /// dynamical variable z
        /// </summary>
        public Double z { get; private set; }

        /// <summary>
        /// initial value of the variable z
        /// </summary>
        public Double z_init { get; private set; }

        /// <summary>
        /// KTzNeuron constructor
        /// </summary>
        /// <param name="K">parameter K</param>
        /// <param name="T">parameter T</param>
        /// <param name="d">parameter delta</param>
        /// <param name="l">parameter lambda</param>
        /// <param name="xR">parameter xR</param>
        /// <param name="x_init">initial value of x</param>
        /// <param name="y_init">initial value of y</param>
        /// <param name="z_init">initial value of z</param>
        public KTzNeuronAvgInp(Int32 ind, Double K, Double T, Double d, Double l, Double xR, Double x_init, Double y_init, Double z_init)
            : base(ind, K, T, x_init, y_init)
        {
            this.d = d;
            this.l = l;
            this.xR = xR;
            this.z_init = this.z = z_init;
        }

        /// <summary>
        /// evaluates this neuron for one timestep
        /// </summary>
        public override void Evolve()
        {
            this.SumInputSignal();
            this.x_prev = this.x;
            this.x = Math.Tanh((this.x - this.K * this.y + this.z + this.Isyn / this.Input.Count) / this.T);
            this.y = this.x_prev;
            this.z = (1.0 - this.d) * this.z - this.l * (this.x_prev - this.xR);
        }

        /// <summary>
        /// evaluates this neuron one timestep with an external stimulus I
        /// </summary>
        /// <param name="I">external stimulus</param>
        public override void Evolve(Double I)
        {
            this.SumInputSignal();
            this.x_prev = this.x;
            this.x = Math.Tanh((this.x - this.K * this.y + this.z + I + this.Isyn / (this.Input.Count + 1)) / this.T);
            this.y = this.x_prev;
            this.z = (1.0 - this.d) * this.z - this.l * (this.x_prev - this.xR);
        }

        /// <summary>
        /// resets neuron variables to the values contained in ic array
        /// </summary>
        /// <param name="ic">new values for neurons variables</param>
        public override void ResetIC(Double[] ic)
        {
            base.ResetIC(ic);
            this.z_init = ic[2];
            this.z = ic[2];
        }
    }


    /// <summary>
    /// implements the neuron of the reference Modeling Neurons by Simple Maps
    /// </summary>
    public class KTLogSimple : KTNeuron
    {
        /// <summary>
        /// parameter H
        /// </summary>
        public Double H { get; private set; }

        /// <summary>
        /// KTLogSimple constructor
        /// </summary>
        /// <param name="K">K parameter</param>
        /// <param name="T">T parameter</param>
        /// <param name="H">H parameter</param>
        /// <param name="x_init">x (membrane potential) initial value</param>
        /// <param name="y_init">y (recurrent variable) initial value</param>
        public KTLogSimple(Int32 ind, Double K, Double T, Double H, Double x_init, Double y_init)
            : base(ind, K, T, x_init, y_init)
        {
            this.H = H;
        }

        private Double logFun(Double u)
        {
            return u / (1 + (u > 0.0D ? u : -u));
        }

        /// <summary>
        /// evaluates this neuron for one timestep
        /// </summary>
        public override void Evolve()
        {
            this.SumInputSignal();
            this.x_prev = this.x;
            this.x = logFun((this.x - this.K * this.y + this.H + this.Isyn) / this.T);
            this.y = this.x_prev;
        }

        /// <summary>
        /// evaluates this neuron one timestep with an external stimulus I
        /// </summary>
        /// <param name="I">external stimulus</param>
        public override void Evolve(Double I)
        {
            this.SumInputSignal();
            this.x_prev = this.x;
            this.x = logFun((this.x - this.K * this.y + this.H + this.Isyn + I) / this.T);
            this.y = this.x_prev;
        }
    }

    /// <summary>
    /// implements the neuron of the reference Modeling Neurons by Simple Maps
    /// </summary>
    public class KTSimple : KTNeuron
    {
        /// <summary>
        /// parameter H
        /// </summary>
        public Double H { get; private set; }

        /// <summary>
        /// KTSimple constructor
        /// </summary>
        /// <param name="K">K parameter</param>
        /// <param name="T">T parameter</param>
        /// <param name="H">H parameter</param>
        /// <param name="x_init">x (membrane potential) initial value</param>
        /// <param name="y_init">y (recurrent variable) initial value</param>
        public KTSimple(Int32 ind, Double K, Double T, Double H, Double x_init, Double y_init)
            : base(ind, K, T, x_init, y_init)
        {
            this.H = H;
        }

        /// <summary>
        /// evaluates this neuron for one timestep
        /// </summary>
        public override void Evolve()
        {
            this.SumInputSignal();
            this.x_prev = this.x;
            this.x = Math.Tanh((this.x - this.K * this.y + this.H + this.Isyn) / this.T);
            this.y = this.x_prev;
        }

        /// <summary>
        /// evaluates this neuron one timestep with an external stimulus I
        /// </summary>
        /// <param name="I">external stimulus</param>
        public override void Evolve(Double I)
        {
            this.SumInputSignal();
            this.x_prev = this.x;
            this.x = Math.Tanh((this.x - this.K * this.y + this.H + this.Isyn + I) / this.T);
            this.y = this.x_prev;
        }
    }

    /// <summary>
    /// abstract class that can be extended to handle one of the following KT Neurons:
    /// MHR Tragtenberg, O Kinouchi: Modeling Neurons by Simple Maps. Journal of Biffurcation and Chaos.
    /// or
    /// SM Kuva, AC Roque, MHR Tragtenberg, O Kinouchi: A minimal model for excitable and bursting elements. Neurocomputing.
    /// </summary>
    public abstract class KTNeuron : INeuron //WithNeighbour<KTSynapse>/*, INeuron*/
    {
        /// <summary>
        /// the index of this neuron (a number that identifies it within a network)
        /// </summary>
        public Int32 Index { get; protected set; }

        /// <summary>
        /// the list of this neuron's neighbours (only the ones that influences this one)
        /// </summary>
        public List<ISynapse> Input { get; private set; }

        /// <summary>
        /// the sum over all neighbours output signals
        /// </summary>
        public Double Isyn { get; private set; }

        /// <summary>
        /// parameter K
        /// </summary>
        public Double K { get; private set; }

        /// <summary>
        /// parameter T
        /// </summary>
        public Double T { get; private set; }

        /// <summary>
        /// previous value of x dynamical variable
        /// </summary>
        protected Double x_prev { get; set; }

        /// <summary>
        /// dynamical variable - membrane potential
        /// </summary>
        public Double x { get; protected set; }

        /// <summary>
        /// dynamical variable - recurrent variable for membrane potential
        /// </summary>
        public Double y { get; protected set; }

        /// <summary>
        /// x initial value
        /// </summary>
        public Double x_init { get; private set; }

        /// <summary>
        /// y initial value
        /// </summary>
        public Double y_init { get; private set; }

        /// <summary>
        /// KTNeuron constructor
        /// </summary>
        /// <param name="K">K parameter</param>
        /// <param name="T">T parameter</param>
        /// <param name="x_init">x (membrane potential) initial value</param>
        /// <param name="y_init">y (recurrent variable) initial value</param>
        protected KTNeuron(Int32 ind, Double K, Double T, Double x_init, Double y_init)
        {
            this.Index = ind;
            this.K = K;
            this.T = T;
            this.x_init = this.x = this.x_prev = x_init;
            this.y_init = this.y = y_init;
            this.Input = new List<ISynapse>();
        }

        /// <summary>
        /// initializes the neighbours with the specified list
        /// </summary>
        /// <param name="s">the list of input synapses</param>
        public void AddInput(List<ISynapse> s)
        {
            this.Input = s;
        }

        /// <summary>
        /// adds a neighbour to this neuron...
        /// </summary>
        /// <param name="s">the input synapse</param>
        public void AddInput(ISynapse s)
        {
            this.Input.Add(s);
        }

        /// <summary>
        /// sums all the neighbours output signals and stores it at this.Isyn
        /// </summary>
        protected void SumInputSignal()
        {
            Int32 i = 0, n = this.Input.Count;
            this.Isyn = 0.0;
            while (i < n)
            {
                this.Isyn += this.Input[i].GetSignal();
                i++;
            }
        }

        /// <summary>
        /// returns the membrane potential of the neuron
        /// </summary>
        /// <returns>this neuron membrane potential</returns>
        public Double GetMembranePotential()
        {
            return this.x;
        }

        /// <summary>
        /// evaluates all this neuron's variables for one timestep
        /// </summary>
        public abstract void Evolve();

        /// <summary>
        /// evaluates this neuron's dynamical variables with external stimulus I
        /// </summary>
        /// <param name="I">the external stimulus on the membrane potential</param>
        public abstract void Evolve(Double I);

        /// <summary>
        /// resets neuron variables to the values contained in ic array
        /// </summary>
        /// <param name="ic">new values for neurons variables</param>
        public virtual void ResetIC(Double[] ic)
        {
            this.x_init = ic[0];
            this.y_init = ic[1];
            this.x = ic[0];
            this.y = ic[1];
        }

        /// <summary>
        /// checks if this neuron is emitting a spike in its current time step iteration and returns true if it is
        /// 
        /// here, a spike is detected whenever the membrane potential crosses x=0 from a negative toa positive value
        /// </summary>
        /// <returns>returns true if the neuron is currently spiking, false otherwise</returns>
        public virtual bool SpikeDetector()
        {
            return (this.x * this.x_prev < 0.0) && (this.x_prev < this.x);
        }
    }

    public class SIElement : ThresholdElement
    {
        public SIElement(Int32 ind, Double Theta, Double x_init)
            : base(ind, Theta, x_init)
        {
        }

        /// <summary>
        /// evaluates all this neuron's variables for one timestep
        /// </summary>
        public override void Evolve()
        {
            this.SumInputSignal();
            if ((this.Isyn / this.TotalInputWeight) >= this.Theta)
            {
                this.x = 1.0;
            }
        }

        /// <summary>
        /// evaluates this neuron's dynamical variables with external stimulus I
        /// </summary>
        /// <param name="I">the external stimulus on the membrane potential</param>
        public override void Evolve(Double I)
        {
            this.SumInputSignal();
            if ((I + (this.Isyn / this.TotalInputWeight)) >= this.Theta)
            {
                this.x = 1.0;
            }
        }
    }

    public abstract class ThresholdElement : INeuron
    {
        /// <summary>
        /// the index of this neuron (a number that identifies it within a network)
        /// </summary>
        public Int32 Index { get; protected set; }

        /// <summary>
        /// the list of this neuron's neighbours (only the ones that influences this one)
        /// </summary>
        public List<ISynapse> Input { get; private set; }

        /// <summary>
        /// the sum over all neighbours output signals
        /// </summary>
        public Double Isyn { get; private set; }

        /// <summary>
        /// the sum over all neighbours connection weight
        /// </summary>
        public Double TotalInputWeight { get; private set; }

        /// <summary>
        /// parameter threshold
        /// </summary>
        public Double Theta { get; private set; }

        /// <summary>
        /// state variable (0 or 1)
        /// </summary>
        public Double x { get; protected set; }

        /// <summary>
        /// ThresholdElement constructor
        /// </summary>
        /// <param name="ind">this neuron index in the network</param>
        /// <param name="Theta">Threshold of this neuron</param>
        /// <param name="x_init">x (membrane potential) initial value</param>
        protected ThresholdElement(Int32 ind, Double Theta, Double x_init)
        {
            this.Index = ind;
            this.Theta = Theta;
            this.x = x_init;
            this.TotalInputWeight = 0.0D;
            this.Input = new List<ISynapse>();
        }

        private void SetTotalInputWeight()
        {
            this.TotalInputWeight = 0.0D;
            Int32 i = 0;
            while (i < this.Input.Count)
            {
                this.TotalInputWeight += this.Input[i].GetCoupling();
                i++;
            }
        }

        /// <summary>
        /// initializes the neighbours with the specified list
        /// </summary>
        /// <param name="s">the list of input synapses</param>
        public void AddInput(List<ISynapse> s)
        {
            this.Input = s;
            this.SetTotalInputWeight();
        }

        /// <summary>
        /// adds a neighbour to this neuron...
        /// </summary>
        /// <param name="s">the input synapse</param>
        public void AddInput(ISynapse s)
        {
            this.Input.Add(s);
            this.SetTotalInputWeight();
        }

        /// <summary>
        /// sums all the neighbours output signals and stores it at this.Isyn
        /// </summary>
        protected void SumInputSignal()
        {
            Int32 i = 0, n = this.Input.Count;
            this.Isyn = 0.0;
            while (i < n)
            {
                this.Isyn += this.Input[i].GetSignal();
                i++;
            }
        }

        /// <summary>
        /// returns the membrane potential of the neuron
        /// </summary>
        /// <returns>this neuron membrane potential</returns>
        public Double GetMembranePotential()
        {
            return this.x;
        }

        /// <summary>
        /// evaluates all this neuron's variables for one timestep
        /// </summary>
        public abstract void Evolve();

        /// <summary>
        /// evaluates this neuron's dynamical variables with external stimulus I
        /// </summary>
        /// <param name="I">the external stimulus on the membrane potential</param>
        public abstract void Evolve(Double I);

        /// <summary>
        /// resets neuron variables to the values contained in ic array
        /// </summary>
        /// <param name="ic">new values for neurons variables</param>
        public virtual void ResetIC(Double[] ic)
        {
            this.x = ic[0];
        }

        /// <summary>
        /// checks if this neuron is emitting a spike in its current time step iteration and returns true if it is
        /// </summary>
        /// <returns>returns true if the neuron is currently spiking, false otherwise</returns>
        public virtual bool SpikeDetector()
        {
            return this.x == 1.0;
        }
    }

    /// <summary>
    /// common interface of the neurons
    /// </summary>
    public interface INeuron
    {
        /// <summary>
        /// the index of this neuron (a number that identifies it within a network)
        /// </summary>
        Int32 Index { get; }

        /// <summary>
        /// main state variable
        /// </summary>
        Double x { get; }

        /// <summary>
        /// the list of this neuron input synapses
        /// </summary>
        List<ISynapse> Input { get; }

        /// <summary>
        /// adds a neighbour to this neuron
        /// </summary>
        /// <param name="s">the input synapse</param>
        void AddInput(ISynapse s);

        /// <summary>
        /// adds a neighbour to this neuron
        /// </summary>
        /// <param name="s">a list of input synapses</param>
        void AddInput(List<ISynapse> s);

        /// <summary>
        /// returns the membrane potential of the neuron
        /// </summary>
        Double GetMembranePotential();

        /// <summary>
        /// evaluates all this neuron's variables for one timestep
        /// </summary>
        void Evolve();

        /// <summary>
        /// evaluates all this neuron's variables with external input I
        /// </summary>
        /// <param name="I">the external input in this neuron's membrane potential</param>
        void Evolve(Double I);

        /// <summary>
        /// resets neuron variables to the values contained in ic array
        /// </summary>
        /// <param name="ic">new values for neurons variables</param>
        void ResetIC(Double[] ic);

        /// <summary>
        /// checks if this neuron is emitting a spike in its current time step iteration and returns true if it is
        /// </summary>
        /// <returns>returns true if the neuron is currently spiking, false otherwise</returns>
        bool SpikeDetector();
    }
}

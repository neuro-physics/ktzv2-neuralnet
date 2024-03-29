using System;
using KTzV2.Neurons;
using KTzV2.Maths.Random;

namespace KTzV2.Synapses
{
    /// <summary>
    /// improvised enum to select synapse
    /// </summary>
    public enum SynapseType
    {
        KTChemicalSynapse,
        KTNoisyChemicalSynapse,
        KTDynamicChemicalSynapse,
        GapJunction,
        PulseCoupling,
        RectifyingGapJunction,
        NormalizedPulseCoupling
    }

    public enum NoiseType
    {
        ProportionalAmplitude,
        GreaterThanJ
    }

    public class SynapseParam
    {
        public INeuron preSynNeuron { get; private set; }
        public INeuron postSynNeuron { get; private set; }
        public Double J { get; private set; }
        public Double J0 { get; private set; }
        public Double noiseAmp { get; private set; }
        public Double noiseRatio { get; private set; }
        public NoiseType noiseType { get; private set; }
        public Double tauf { get; private set; }
        public Double taug { get; private set; }
        public UInt64 seed { get; private set; }
        public Double alpha { get; private set; }
        public Double u { get; private set; }
        public Double tauJ { get; private set; }
        public Double dt { get; private set; }
        public SynapseParam(INeuron preNeu, INeuron postNeu, Double tauf, Double taug, Double J, Double noiseAmp, Double noiseRatio, UInt64 seed, NoiseType noiseType, Double J0, Double tauJ, Double alpha, Double u, Double dt)
        {
            this.preSynNeuron = preNeu;
            this.postSynNeuron = postNeu;
            this.tauf = tauf;
            this.taug = taug;
            this.J = J;
            this.noiseAmp = noiseAmp;
            this.noiseRatio = noiseRatio;
            this.seed = seed;
            this.noiseType = noiseType;
            this.J0 = J0;
            this.tauJ = tauJ;
            this.alpha = alpha;
            this.u = u;
            this.dt = dt;
        }
    }

    public static class SynapseFactory
    {
        public static ISynapse GetSynapse(SynapseType synType, SynapseParam par)
        {
            switch (synType)
            {
                case SynapseType.KTChemicalSynapse:
                    return new KTChemicalSynapse(par.preSynNeuron, par.postSynNeuron, par.J, par.tauf, par.taug);
                case SynapseType.KTNoisyChemicalSynapse:
                    return new KTNoisyChemicalSynapse(par.preSynNeuron, par.postSynNeuron, par.J, par.tauf, par.taug, par.noiseAmp, par.noiseRatio, par.seed, par.noiseType);
                case SynapseType.RectifyingGapJunction:
                case SynapseType.GapJunction:
                    return new GapJunction(par.preSynNeuron, par.postSynNeuron, par.J, synType == SynapseType.RectifyingGapJunction);
                case SynapseType.KTDynamicChemicalSynapse:
                    return new KTDynamicChemicalSynapse(par.preSynNeuron, par.postSynNeuron, par.J0, par.tauf, par.taug, par.alpha, par.u, par.tauJ, par.dt);
                case SynapseType.NormalizedPulseCoupling:
                case SynapseType.PulseCoupling:
                    return new PulseCoupling(par.preSynNeuron, par.J, synType == SynapseType.NormalizedPulseCoupling);
                default:
                    throw new ArgumentOutOfRangeException("Unrecognized synapse type! " + synType.ToString());
            }
        }
    }

    /// <summary>
    /// class to handle an output synapse according to the reference:
    /// SM Kuva, AC Roque, MHR Tragtenberg, O Kinouchi: A minimal model for excitable and bursting elements. Neurocomputing.
    /// </summary>
    public class KTChemicalSynapse : KTSynapse
    {
        /// <summary>
        /// dynamical variable with output signal
        /// </summary>
        public Double f { get; private set; }

        /// <summary>
        /// dynamical auxiliary variable
        /// </summary>
        private Double g { get; set; }

        /// <summary>
        /// parameter - time constant for f variable
        /// </summary>
        public Double oneOverTau_f { get; private set; }

        /// <summary>
        /// parameter - time constant for g variable
        /// </summary>
        public Double oneOverTau_g { get; private set; }

        /// <summary>
        /// step function to decide if this synapse should act (for membrane potential of the presyn neuron > 0)
        /// </summary>
        private Double thetaJ { get; set; }

        public KTChemicalSynapse(INeuron nPre, INeuron nPost, Double J, Double tau_f, Double tau_g)
            : base(nPre, nPost, J)
        {
            this.oneOverTau_f = 1.0 / tau_f;
            this.oneOverTau_g = 1.0 / tau_g;
            this.f = 0.0;
            this.g = 0.0;
        }

        /// <summary>
        /// evaluates this synaptic signal for one timestep
        /// </summary>
        public override void EvolveLocal()
        {
            this.thetaJ = (this.PreSynapticNeuron.GetMembranePotential() > 0.0 ? this.J : 0.0);
            this.f = (1.0 - this.oneOverTau_f) * this.f + this.g;
            this.g = (1.0 - this.oneOverTau_g) * this.g + this.thetaJ;
        }

        /// <summary>
        /// gets the signal that this synapse generates (the f variable)
        /// </summary>
        /// <returns>f variable - a double value with the signal of this synapse</returns>
        public override Double GetSignal()
        {
            return this.f;
        }

        /// <summary>
        /// Resets output synaptic signal
        /// </summary>
        public override void ResetSignal()
        {
            this.f = 0.0D;
            this.g = 0.0D;
        }
    }

    /// <summary>
    /// class to handle an output synapse according to the reference:
    /// SM Kuva, AC Roque, MHR Tragtenberg, O Kinouchi: A minimal model for excitable and bursting elements. Neurocomputing.
    /// </summary>
    public class KTDynamicChemicalSynapse : KTSynapse
    {
        /// <summary>
        /// dynamical variable with output signal
        /// </summary>
        public Double f { get; private set; }

        /// <summary>
        /// dynamical auxiliary variable
        /// </summary>
        private Double g { get; set; }

        /// <summary>
        /// parameter - time constant for f variable
        /// </summary>
        public Double oneOverTau_f { get; private set; }

        /// <summary>
        /// parameter - time constant for g variable
        /// </summary>
        public Double oneOverTau_g { get; private set; }

        /// <summary>
        /// parameter - time constant of the coupling parameter
        /// </summary>
        public Double oneOverTau_J { get; private set; }

        /// <summary>
        /// parameter - maximum amplitude of the synapse
        /// </summary>
        public Double alphaOverU { get; private set; }

        /// <summary>
        /// parameter - amount of coupling lost when a spike occurs
        /// </summary>
        public Double u { get; private set; }

        /// <summary>
        /// step function to decide if this synapse should act (for membrane potential of the presyn neuron > 0)
        /// </summary>
        private Double thetaJ { get; set; }

        /// <summary>
        /// Creates a KT Dynamic Chemical Synapse
        /// </summary>
        /// <param name="nPre">the pre-syn neuron</param>
        /// <param name="J0">the initial value of the coupling</param>
        /// <param name="tau_f">the time constant of f variable</param>
        /// <param name="tau_g">the time constant of g variable</param>
        /// <param name="alpha">the alpha parameter: alpha/u is the maximum value of J (the coupling parameter)</param>
        /// <param name="u">the u parameter: the decrease on the coupling when a spike occurs</param>
        /// <param name="tau_J">the time constant of J (coupling). It should be of the order of tau*N, where N is the number of neurons and tau is the time constant of the stimulus</param>
        /// <param name="dt">the timestep to rescale the parameters, as this is a discretization of a differential equation</param>
        public KTDynamicChemicalSynapse(INeuron nPre, INeuron nPost, Double J0, Double tau_f, Double tau_g, Double alpha, Double u, Double tau_J, Double dt)
            : base(nPre, nPost, J0)
        {
            this.oneOverTau_f = 1.0 / tau_f;
            this.oneOverTau_g = 1.0 / tau_g;
            this.oneOverTau_J = dt / tau_J; // rescaling this term according to dt
            this.alphaOverU = alpha / u; // this term needs not to be rescaled because it is a constant which is multiplied by dt when it is multiplied by oneOverTau_J
            this.u = dt * u; // rescaling this term according to dt
            this.f = 0.0;
            this.g = 0.0;
        }

        /// <summary>
        /// evaluates this synaptic signal for one timestep
        /// </summary>
        public override void EvolveLocal()
        {
            this.thetaJ = this.J + (this.PreSynapticNeuron.GetMembranePotential() > 0.0 ? (this.alphaOverU - this.J) * this.oneOverTau_J : (this.alphaOverU - this.J) * this.oneOverTau_J - this.u * this.J);
            this.f = (1.0 - this.oneOverTau_f) * this.f + this.g;
            this.g = (1.0 - this.oneOverTau_g) * this.g + this.thetaJ;
        }

        /// <summary>
        /// gets the signal that this synapse generates (the f variable)
        /// </summary>
        /// <returns>f variable - a double value with the signal of this synapse</returns>
        public override Double GetSignal()
        {
            return this.f;
        }

        /// <summary>
        /// Resets output synaptic signal
        /// </summary>
        public override void ResetSignal()
        {
            this.f = 0.0D;
            this.g = 0.0D;
        }
    }

    /// <summary>
    /// class to handle an output chemical synapse with noise. The chemical synapse is according to the reference:
    /// SM Kuva, AC Roque, MHR Tragtenberg, O Kinouchi: A minimal model for excitable and bursting elements. Neurocomputing.
    /// </summary>
    public class KTNoisyChemicalSynapse : KTSynapse
    {
        /// <summary>
        /// dynamical variable with output signal
        /// </summary>
        public Double f { get; private set; }

        /// <summary>
        /// dynamical auxiliary variable
        /// </summary>
        private Double g { get; set; }

        /// <summary>
        /// parameter - time constant for f variable
        /// </summary>
        public Double oneOverTau_f { get; private set; }

        /// <summary>
        /// parameter - time constant for g variable
        /// </summary>
        public Double oneOverTau_g { get; private set; }

        /// <summary>
        /// step function to decide if this synapse should act (for membrane potential of the presyn neuron > 0)
        /// </summary>
        private Double thetaJ { get; set; }

        /// <summary>
        /// random number to generate noise over J
        /// </summary>
        private HomogeneousRand Jnoise { get; set; }

        private Func<Double> GetJ { get; set; }

        //private Int32 JSign { get; set; }

        /// <summary>
        /// the bounds of the noise signal
        /// </summary>
        private Double NoiseAmp { get; set; }

        public KTNoisyChemicalSynapse(INeuron nPre, INeuron nPost, Double J, Double tau_f, Double tau_g, Double noiseAmp, Double noiseRatio, UInt64 seed, NoiseType nt)
            : base(nPre, nPost, J)
        {
            this.Jnoise = new HomogeneousRand(seed);

            if (nt == NoiseType.ProportionalAmplitude)
            {
                this.NoiseAmp = Math.Abs(J * noiseRatio);
                this.GetJ = this.GetJProportionalAmplitude;
            }
            else if (nt == NoiseType.GreaterThanJ)
            {
                this.NoiseAmp = Math.Abs(noiseAmp) * Math.Sign(J);
                this.GetJ = this.GetJGreaterThanJ;
            }
            else
            {
                throw new ArgumentOutOfRangeException("Unrecognized NoiseType! " + nt.ToString());
            }

            this.oneOverTau_f = 1.0 / tau_f;
            this.oneOverTau_g = 1.0 / tau_g;
            this.f = 0.0;
            this.g = 0.0;
        }

        private Double GetJProportionalAmplitude()
        {
            return this.J + (this.NoiseAmp * (2.0D * this.Jnoise.GetRandomFull() - 1.0D));
        }

        private Double GetJGreaterThanJ()
        {
            return this.J + (this.NoiseAmp * this.Jnoise.GetRandomFull());
        }

        /// <summary>
        /// evaluates this synaptic signal for one timestep
        /// </summary>
        public override void EvolveLocal()
        {
            // otimizar esta parte... ajustar no construtor this.NoiseAmp = Math.Abs(noiseAmp) * this.JSign
            //this.thetaJ = (this.PreSynapticNeuron.GetMembranePotential() > 0.0 ? this.J + (this.JSign * this.NoiseAmp * this.Jnoise.GetRandomFull()) : 0.0);
            this.thetaJ = (this.PreSynapticNeuron.GetMembranePotential() > 0.0D ? this.GetJ() : 0.0D);
            this.f = (1.0 - this.oneOverTau_f) * this.f + this.g;
            this.g = (1.0 - this.oneOverTau_g) * this.g + this.thetaJ;
        }

        /// <summary>
        /// gets the signal that this synapse generates (the f variable)
        /// </summary>
        /// <returns>f variable - a double value with the signal of this synapse</returns>
        public override Double GetSignal()
        {
            return this.f;
        }

        /// <summary>
        /// Resets output synaptic signal
        /// </summary>
        public override void ResetSignal()
        {
            this.f = 0.0D;
            this.g = 0.0D;
        }
    }

    /// <summary>
    /// common class for all KT synapses
    /// </summary>
    public abstract class KTSynapse : ISynapse//FromPreSynaptic//<INeuron>/*, ISynapse*/
    {
        /// <summary>
        /// the neuron which generates this synapse
        /// </summary>
        public INeuron PreSynapticNeuron { get; protected set; }

        /// <summary>
        /// the neuron which generates this synapse
        /// </summary>
        public INeuron PostSynapticNeuron { get; protected set; }

        /// <summary>
        /// parameter or variable - influence of the presyn neuron over its' neighbours
        /// </summary>
        public Double J { get; protected set; }

        public Action Evolve { get; private set; }

        /// <summary>
        /// constructor of a KTSynapse
        /// </summary>
        /// <param name="nPre">the presynaptic neuron</param>
        /// <param name="nPost">the postsynaptic neuron</param>
        /// <param name="J">the coupling of the presynaptic neuron with its neighbours</param>
        protected KTSynapse(INeuron nPre, INeuron nPost, Double J)
        {
            this.J = J;
            this.PreSynapticNeuron = nPre;
            this.PostSynapticNeuron = nPost;
            this.Evolve = this.EvolveLocal;
        }

        /// <summary>
        /// evaluating method - evaluates one timestep of this synapse
        /// </summary>
        public abstract void EvolveLocal();

        /// <summary>
        /// gets the signal that this synapse generates
        /// </summary>
        /// <returns>a double value with the signal of this synapse</returns>
        public abstract Double GetSignal();

        public Double GetCoupling()
        {
            return this.J;
        }

        /// <summary>
        /// Resets output synaptic signal
        /// </summary>
        public abstract void ResetSignal();
    }

    /// <summary>
    /// creates a Gap Junction synapse
    /// according to:
    /// Deterministic excitable media under Poisson drive: Power law responses,
    /// spiral waves, and dynamic range
    /// PHYSICAL REVIEW E77, 051911 2008
    /// Tiago L. Ribeiro and Mauro Copelli
    /// </summary>
    public class GapJunction : ISynapse
    {
        /// <summary>
        /// the presyn neuron
        /// </summary>
        public INeuron PreSynapticNeuron { get; private set; }

        /// <summary>
        /// the post syn neuron
        /// </summary>
        public INeuron PostSynapticNeuron { get; private set; }

        /// <summary>
        /// channel conductance
        /// </summary>
        public Double J { get; private set; }

        /// <summary>
        /// stores the state of the synapse according to neuron states in the previous timestep in order to correct interaction timing
        /// </summary>
        private Double f { get; set; }

        public Action Evolve { get; private set; }

        private Double dV { get; set; }

        /// <summary>
        /// constructor of a KTSynapse
        /// </summary>
        /// <param name="nPre">the presynaptic neuron</param>
        /// <param name="nPost">the postsynaptic neuron</param>
        /// <param name="J">the coupling (conductance of the channel) of the presynaptic neuron with its neighbours</param>
        /// <param name="rectifying">if true, then only activates when the presynaptic potential is greater than postsynaptic potential</param>
        public GapJunction(INeuron nPre, INeuron nPost, Double J, bool rectifying=false)
        {
            this.J = J;
            this.PreSynapticNeuron = nPre;
            this.PostSynapticNeuron = nPost;
            this.f = 0.0;
            if (rectifying)
            {
                this.Evolve = this.EvolveRectifying;
            }
            else
            {
                this.Evolve = this.EvolveNonRectifying;
            }
        }

        /// <summary>
        /// evolve a rectifying gap junction -- Erik de Schutter book (Calabrese & Prinz, p.288)
        /// </summary>
        public void EvolveRectifying()
        {
            // I'm doing pre - post because of the "-" sign in the definition of Calabrese & Prinz, eq 12.1, book Computational Modeling for Neuroscientists, p. 286
            // in this way, it must enter as a "+" term in the sum of currents of the "dV/dT" equation
            // for the same reason, I inverted the sign to be < 0.0 in the condition check (instead of > 0.0)
            dV     = this.PreSynapticNeuron.GetMembranePotential() - this.PostSynapticNeuron.GetMembranePotential();
            this.f = this.J * dV * (dV < 0.0D ? 0.0D : 1.0D);
        }

        /// <summary>
        /// evolve a non-rectifying gap junction -- Erik de Schutter book (Calabrese & Prinz, p.288)
        /// </summary>
        public void EvolveNonRectifying()
        {
            // I'm doing pre - post because of the "-" sign in the definition of Calabrese & Prinz, eq 12.1, book Computational Modeling for Neuroscientists, p. 286
            // in this way, it must enter as a "+" term in the sum of currents of the "dV/dT" equation
            this.f = this.J * (this.PreSynapticNeuron.GetMembranePotential() - this.PostSynapticNeuron.GetMembranePotential());
        }

        /// <summary>
        /// gets synaptic signal j->i
        /// </summary>
        /// <returns>J * (Vj - Vi)</returns>
        public Double GetSignal()
        {
            return this.f;
        }

        public Double GetCoupling()
        {
            return this.J;
        }

        /// <summary>
        /// Resets output synaptic signal
        /// </summary>
        public void ResetSignal()
        {
            return;
        }
    }

    /// <summary>
    /// pulse coupling
    /// I(t) = J*V_pre
    /// </summary>
    public class PulseCoupling : ISynapse
    {
        /// <summary>
        /// a reference to the neuron which generates this synapse
        /// </summary>
        public INeuron PreSynapticNeuron { get; private set; }

        /// <summary>
        /// a reference to the neuron which reads this signal
        /// </summary>
        public INeuron PostSynapticNeuron { get; private set; }

        /// <summary>
        /// parameter - influence of the presyn neuron over its' neighbours
        /// </summary>
        public Double J { get; private set; }

        /// <summary>
        /// stores the state of the synapse according to neuron states in the previous timestep in order to correct interaction timing
        /// </summary>
        private Double f { get; set; }

        public Action Evolve { get; private set; }

        /// <summary>
        /// pulse coupling
        /// </summary>
        /// <param name="nPre">presynaptic neuron</param>
        /// <param name="J">conductance of the channel</param>
        public PulseCoupling(INeuron nPre, Double J, bool normalizePulse = false)
        {
            this.PreSynapticNeuron = nPre;
            this.J = J;
            this.f = 0.0;
            if (normalizePulse)
            {
                this.Evolve = this.EvolveLocalNormalized;
            }
            else
            {
                this.Evolve = this.EvolveLocal;
            }
        }

        /// <summary>
        /// evaluates one timestep of this synapse
        /// </summary>
        public void EvolveLocalNormalized()
        {
            this.f = this.J * (this.PreSynapticNeuron.GetMembranePotential()>0?1.0D:0.0D);
        }

        public void EvolveLocal()
        {
            this.f = this.J * this.PreSynapticNeuron.GetMembranePotential();
        }

        /// <summary>
        /// gets the signal that this synapse generates
        /// </summary>
        /// <returns>a double value with the signal of this synapse</returns>
        public Double GetSignal()
        {
            return this.f;
        }

        public Double GetCoupling()
        {
            return this.J;
        }

        /// <summary>
        /// Resets output synaptic signal
        /// </summary>
        public void ResetSignal()
        {
            return;
        }
    }

    public interface ISynapse
    {
        /// <summary>
        /// a reference to the neuron which generates this synapse
        /// </summary>
        INeuron PreSynapticNeuron { get; }

        /// <summary>
        /// a reference to the neuron which reads this signal
        /// </summary>
        INeuron PostSynapticNeuron { get; }

        /// <summary>
        /// parameter - influence of the presyn neuron over its' neighbours
        /// </summary>
        Double J { get; }

        /// <summary>
        /// evaluates one timestep of this synapse
        /// </summary>
        Action Evolve { get; }

        /// <summary>
        /// gets the signal that this synapse generates
        /// </summary>
        /// <returns>a double value with the signal of this synapse</returns>
        Double GetSignal();

        /// <summary>
        /// gets the coupling intensity (J) of this synapse
        /// </summary>
        /// <returns>value of J parameter</returns>
        Double GetCoupling();

        /// <summary>
        /// Resets output synaptic signal
        /// </summary>
        void ResetSignal();
    }
}

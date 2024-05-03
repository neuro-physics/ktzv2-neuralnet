using System;
using KTzV2.Maths.Random;

namespace KTzV2.Stimuli
{
    public enum StimulusType
    {
        Delta,
        DeltaTrain,
        PoissonProcess,
        DeltaWhenInactive,
        None
    }

    public static class FExternalStimulus
    {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="t0">time for first stimulus</param>
        /// <param name="stimulusIntensity">the (mean) intensity of the stimulus</param>
        /// <param name="deltat">time interval between two consecutive delta stimulus</param>
        /// <param name="n">amount of delta stimulus</param>
        /// <param name="r">the inverse of the Poisson rate (timesteps/event)</param>
        /// <param name="intensityStdDev">the std deviation of the generated signal</param>
        /// <returns>The desired stimulus type</returns>
        public static ExternalStimulusBase Get(StimulusType st, Int32 t0, Double stimulusIntensity, Int32 deltat, Int32 n, Double r, Double intensityStdDev)
        {
            switch (st)
            {
                case StimulusType.Delta: return new DeltaStimulus(t0, stimulusIntensity);
                case StimulusType.DeltaTrain: return new DeltaTrainStimulus(t0, stimulusIntensity, deltat, n);
                case StimulusType.PoissonProcess: return new PoissonProcess(r, stimulusIntensity, intensityStdDev);
                case StimulusType.DeltaWhenInactive: return new DeltaWhenInactive(stimulusIntensity);
                case StimulusType.None: return new ExternalStimulusBase();
            }
            throw new ArgumentOutOfRangeException("unrecognized stimulus type!");
        }
    }

    /// <summary>
    /// interface to implement an external stimulus
    /// </summary>
    public interface IExternalStimulus
    {
        /// <summary>
        /// timestep for stimulation
        /// </summary>
        Int32 timeForStimulus { get; }

        /// <summary>
        /// gets the signal of the stimulus for the given timestep
        /// </summary>
        /// <param name="timestep">the timestep to get the signal</param>
        /// <returns></returns>
        Double GetSignal(Int32 timestep);
    }

    public class ExternalStimulusBase : IExternalStimulus
    {
        /// <summary>
        /// timestep for stimulation
        /// </summary>
        public Int32 timeForStimulus { get; protected set; }

        public ExternalStimulusBase()
        { }
        
        /// <summary>
        /// gets the signal of the stimulus for the given timestep
        /// </summary>
        /// <param name="timestep">the timestep to get the signal</param>
        /// <returns>the intensity of the signal of the stimulus</returns>
        public virtual Double GetSignal(Int32 timestep)
        {
            return 0.0;
        }

        public virtual void SetNextStimulusTime(Int32 timeForNextStimulus)
        {
            return;
        }
    }

    public class DeltaWhenInactive : ExternalStimulusBase
    {
        private Double StimulusIntensity { get; set; }
        //public Int32 timeForStimulus { get; private set; }

        public DeltaWhenInactive(Double stimulusIntensity)
        {
            this.StimulusIntensity = stimulusIntensity;
        }

        public override Double GetSignal(Int32 timestep)
        {
            if (timestep == this.timeForStimulus)
            {
                return this.StimulusIntensity;
            }
            else
            {
                return 0.0;
            }
        }

        public override void SetNextStimulusTime(Int32 timeForNextStimulus)
        {
            this.timeForStimulus = timeForNextStimulus;
        }
    }

    public class DeltaTrainStimulus : ExternalStimulusBase
    {
        private Double stimulusIntensity { get; set; }

        private Int32[] timeForStimulusDT { get; set; }

        private Int32 nextTime { get; set; }

        /// <summary>
        /// creates a delta train stimulus (a set of n delta stimulus, each separated by deltat time of another)
        /// </summary>
        /// <param name="t0">time for first stimulus</param>
        /// <param name="deltat">time interval between two consecutive delta stimulus</param>
        /// <param name="n">amount of delta stimulus</param>
        /// <param name="stimulusIntensity">the intensity of the stimulus</param>
        public DeltaTrainStimulus(Int32 t0, Double stimulusIntensity, Int32 deltat, Int32 n)
        {
            if (n < 1) throw new ArgumentOutOfRangeException("n must be >= 1");
            this.stimulusIntensity = stimulusIntensity;
            this.timeForStimulusDT = new Int32[n + 1];
            this.timeForStimulusDT[0] = t0;
            Int32 i = 1;
            while (i < n)
            {
                this.timeForStimulusDT[i] = t0 + i * deltat;
                i++;
            }
            //
            // if this.nextTime == n, there won't be anything in this.timeForStimulus, so the program will throw an exception: ArgumentOutOfRangeException
            // to fix this, we put t0 as this.timeForStimulus[n], so it will never be equal to any timestep after this.timeForStimulus[n-1] and
            // this.nextTime won't be summed up anymore.
            //
            this.timeForStimulusDT[n] = t0;
            this.nextTime = 0;
        }

        public override Double GetSignal(Int32 timestep)
        {
            if (timestep == this.timeForStimulusDT[this.nextTime])
            {
                this.nextTime++;
                return stimulusIntensity;
            }
            else
            {
                return 0.0;
            }
        }
    }

    public class DeltaStimulus : ExternalStimulusBase
    {
        private Double stimulusIntensity { get; set; }

        //private Int32 timeForStimulus { get; set; }

        public DeltaStimulus(Int32 timeForStimulus, Double stimulusIntensity)
        {
            this.stimulusIntensity = stimulusIntensity;
            this.timeForStimulus = timeForStimulus;
        }

        public override Double GetSignal(Int32 timestep)
        {
            if (timestep == this.timeForStimulus)
            {
                return stimulusIntensity;
            }
            else
            {
                return 0.0;
            }
        }
    }

    /// <summary>
    /// Class to create a list of signals (randomly distributed like a gaussian) organized in time like a Poisson process (exponential interval distribution between events)
    /// </summary>
    public class PoissonProcess : ExternalStimulusBase
    {
        /// <summary>
        /// gaussian random number generator to get the signal intensity
        /// </summary>
        private GaussianRand signalIntensity { get; set; }

        /// <summary>
        /// random number to calculate the interval distribution
        /// </summary>
        private HomogeneousRand rand { get; set; }

        /// <summary>
        /// inverse of the Poisson process rate (unit: timesteps/events)
        /// </summary>
        public Double r { get; private set; }

        /// <summary>
        /// creates a Poisson process of undefined time duration
        /// </summary>
        /// <param name="r">the inverse of the Poisson rate (timesteps/event)</param>
        /// <param name="signalMean">the mean intensity of the generated signal</param>
        /// <param name="signalStdDev">the std deviation of the generated signal</param>
        public PoissonProcess(Double r, Double signalMean, Double signalStdDev)
        {
            this.r = 1.0D - Math.Exp(-r); // this is the probability of firing, and it is constant
            this.rand = new HomogeneousRand();
            this.signalIntensity = new GaussianRand(signalMean, signalStdDev);
        }

        /// <summary>
        /// gets the signal intensity at the specified time
        /// </summary>
        /// <param name="timestep">the specified time to get the signal</param>
        /// <returns>the signal intensity (gaussian distributed around its mean)</returns>
        public override Double GetSignal(Int32 timestep)
        {
            return this.signalIntensity.getRandom()*(this.rand.GetRandomGT0LT1() < this.r? 1.0D : 0.0D);
        }
    }

    /*
    /// <summary>
    /// Class to create a list of signals (randomly distributed like a gaussian) organized in time like a Poisson process (exponential interval distribution between events)
    /// </summary>
    /*
    public class PoissonProcess : ExternalStimulusBase
    {
        /// <summary>
        /// gaussian random number generator to get the signal intensity
        /// </summary>
        private GaussianRand signalIntensity { get; set; }

        /// <summary>
        /// random number to calculate the interval distribution
        /// </summary>
        private HomogeneousRand rand { get; set; }

        /// <summary>
        /// inverse of the Poisson process rate (unit: timesteps/events)
        /// </summary>
        public Double r { get; private set; }

        /// <summary>
        /// indicates the time when will occur the next fire of this event
        /// </summary>
        private Int32 timeForNextFire { get; set; }

        ///// <summary>
        ///// the state of the Poisson signal... each index is for each timestep the signal is running
        ///// </summary>
        //private Int32[] signalState { get; set; }

        ///// <summary>
        ///// creates a Poisson process with the specified parameters
        ///// </summary>
        ///// <param name="r">the inverse of the Poisson rate (timesteps/event)</param>
        ///// <param name="signalMean">the mean intensity of the signal to generate</param>
        ///// <param name="signalStdDev">the std dev of the signal to generate</param>
        ///// <param name="totalTimesteps">total amount of timesteps that this process will take</param>
        //public PoissonProcess(Double r, Double signalMean, Double signalStdDev, Int32 totalTimesteps)
        //    : this(r, signalMean, signalStdDev)
        //{
        //    this.GenerateSignal(totalTimesteps);
        //}

        /// <summary>
        /// creates a Poisson process of undefined time duration
        /// </summary>
        /// <param name="r">the inverse of the Poisson rate (timesteps/event)</param>
        /// <param name="signalMean">the mean intensity of the generated signal</param>
        /// <param name="signalStdDev">the std deviation of the generated signal</param>
        public PoissonProcess(Double r, Double signalMean, Double signalStdDev)
        {
            this.r = r;
            this.rand = new HomogeneousRand();
            this.signalIntensity = new GaussianRand(signalMean, signalStdDev);
            this.timeForNextFire = this.NextTimeInterval();
        }

        ///// <summary>
        ///// generates an array of size totalTimesteps with the Poisson signal
        ///// </summary>
        ///// <param name="totalTimesteps">total amount of timesteps that this process will take</param>
        //private void GenerateSignal(Int32 totalTimesteps)
        //{
        //    this.signalState = new Int32[totalTimesteps];
        //    Int32 nextFire = this.NextTimeInterval(), i = 0;
        //    while (i < totalTimesteps)
        //    {
        //        if (i == nextFire)
        //        {
        //            this.signalState[i] = 1;
        //            nextFire += this.NextTimeInterval();
        //        }
        //        else
        //        {
        //            this.signalState[i] = 0;
        //        }
        //        i++;
        //    }
        //}

        /// <summary>
        /// get the next time interval (in timesteps) between two adjacent signals exponentialy distributed
        /// </summary>
        /// <returns>the time interval in timesteps</returns>
        private Int32 NextTimeInterval()
        {
            return 1 - (Int32)(Math.Log(this.rand.GetRandomGT0LT1()) * this.r);
        }

        /// <summary>
        /// gets the signal intensity at the specified time
        /// </summary>
        /// <param name="timestep">the specified time to get the signal</param>
        /// <returns>the signal intensity (gaussian distributed around its mean)</returns>
        public override Double GetSignal(Int32 timestep)
        {
            if (timestep == this.timeForNextFire)
            {
                this.timeForNextFire += this.NextTimeInterval();
                return this.signalIntensity.getRandom();
            }
            else
            {
                return 0.0;
            }
            //return this.signalState[timestep] * this.gRand.getRandom();
        }
    }
    /**/
}

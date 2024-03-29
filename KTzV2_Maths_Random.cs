using System;
using MersenneTwister;

namespace KTzV2.Maths.Random
{
    public class GaussianRand
    {
        private MT19937 rand { get; set; } // the random number
        public Double stdDev { get; private set; }
        public Double mean { get; private set; }
        public Func<Double> getRandom { get; private set; }

        public GaussianRand(Double mean, Double stdDev)
            : this (mean, stdDev, Convert.ToUInt64(DateTime.Now.Ticks)) // seeding the random generator with the current time
        {
            //this.mean = mean;
            //this.stdDev = stdDev;
            //rand = new MT19937();
            //rand.init_genrand(Convert.ToUInt64(DateTime.Now.Ticks)); // seeding the random generator with the current time
        }

        public GaussianRand(Double mean, Double stdDev, UInt64 seed)
        {
            this.mean = mean;
            this.stdDev = stdDev;
            rand = new MT19937();
            rand.init_genrand(seed);
            if (stdDev == 0.0D)
            {
                this.getRandom = this.getRandom0;
            }
            else
            {
                this.getRandom = this.getRandomG;
            }
        }

        private Double getRandom0()
        {
            return this.mean;
        }

        private Double getRandomG()
        {
            Double v1;
            Double v2;
            Double R;

            do
            {
                v1 = 2.0 * rand.genrand_real3() - 1.0;
                v2 = 2.0 * rand.genrand_real3() - 1.0;
                R = v1 * v1 + v2 * v2;
            } while (R >= 1.0);

            return this.mean + this.stdDev * v2 * Math.Sqrt(-2.0 * Math.Log(R) / R);
        }
    }

    public class HomogeneousRand
    {
        private MT19937 rand;
        //private Random rand;
        //private Int32 seed;

        public HomogeneousRand()
            : this (Convert.ToUInt64(DateTime.Now.Ticks))
        {
            //seed = DateTime.Now.Hour + DateTime.Now.Minute + DateTime.Now.Second + DateTime.Now.Millisecond;
            //rand = new MT19937();
            //rand.init_genrand(Convert.ToUInt64(DateTime.Now.Ticks)); // seeding the random generator with the current time
        }

        public HomogeneousRand(UInt64 seed)
        {
            rand = new MT19937();
            rand.init_genrand(seed); // seeding the random generator with the current time
        }

        public Double GetRandomFull()
        {
            //rand = new Random(seed++);
            return rand.genrand_real1();
        }

        public Double GetRandomLT1()
        {
            return rand.genrand_real2();
        }

        public Double GetRandomGT0LT1()
        {
            return rand.genrand_real3();
        }
    }
}
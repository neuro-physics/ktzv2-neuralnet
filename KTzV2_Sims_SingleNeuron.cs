using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using KTzV2.Data.Header;
using KTzV2.Neurons;

namespace KTzV2.Sims.SingleNeuron
{
    class KTzNeuronSimulator
    {
        //KTzNeuron neuron;

        Int32 t0;
        Int32 tTotal;

        List<Double> xData;
        List<Int32> tData;

        public KTzNeuronSimulator()
        {
            this.CreateNeuron();
            this.GetParameters();
        }

        public void Run()
        {
            Int32 t = 0;
            while (t < tTotal)
            {
                //neuron.Evolve();
                t++;
            }
        }

        private void GetParameters()
        {
            this.t0 = KTzHeader.GetPar_Int32(KTzParameters.nStart);
            this.tTotal = KTzHeader.GetPar_Int32(KTzParameters.nSteps) + 1;

            this.xData = new List<Double>();
            this.tData = new List<Int32>();

            System.IO.MemoryStream ms = new System.IO.MemoryStream();
            System.IO.StreamWriter sw = new System.IO.StreamWriter(ms);
            //sw.Write
        }

        private void CreateNeuron()
        {
            //this.neuron = new KTzNeuron(0, KTzHeader.GetPar_Double(KTzParameters.K), KTzHeader.GetPar_Double(KTzParameters.T), KTzHeader.GetPar_Double(KTzParameters.d), KTzHeader.GetPar_Double(KTzParameters.l), KTzHeader.GetPar_Double(KTzParameters.xR), KTzHeader.GetPar_Double(KTzParameters.x0), KTzHeader.GetPar_Double(KTzParameters.y0), KTzHeader.GetPar_Double(KTzParameters.z0));
        }

        /// <summary>
        /// finds the t-interval between two consecutive peaks (spikes) of the function defined by the data t,x[nIndex],
        /// assuming the peak is centered between two consecutive crossings of the x=0 axis.
        /// </summary>
        /// <param name="nIndex">the index of the neuron to find the ISI</param>
        private List<Int32> findISI(Int32 nIndex)
        {
            List<Int32> isiData = new List<Int32>();
            Int32 n = this.xData.Count;
            Int32 i, iplus, crossCounter = 0;
            Int32 t1 = 0, t2 = 0, ts = 0, tsA = 0;

            n--;
            i = 0;
            while (this.xData[i] >= 0.0) // this is needed to start with data that are less than zero
            {
                i++;
            }
            while (i < n)
            {
                iplus = i + 1;
                if (this.xData[i] * this.xData[iplus] < 0)
                {
                    if (crossCounter % 2 == 0)
                    {
                        t1 = (tData[i] + tData[iplus]) / 2;
                    }
                    else
                    {
                        t2 = (tData[i] + tData[iplus]) / 2;
                        ts = (t1 + t2) / 2;
                        isiData.Add(ts - tsA);
                        tsA = ts;
                    }
                    crossCounter++;
                }
                i = iplus;
            }
            if (isiData.Count > 0)
            {
                isiData.RemoveAt(0); // as we start with xs = 0, the first distance calculated is now an ISI distance, but it shouldn't, so we must delete it
            }
            return isiData;
        }

        /// <summary>
        /// counts the peaks of the neuron k between t1 and t2
        /// </summary>
        /// <param name="k">the neuron index</param>
        /// <param name="t1">the first index of time</param>
        /// <param name="t2">the second index of time. if t2 == t1, then we count from t1 to the end of data</param>
        /// <returns>sum over k-th neuron spikes between t1 and t2</returns>
        private Int32 CountNonConsecutiveDataPeaks(Int32 k, Int32 t1, Int32 t2)
        {
            Int32 peakCounter = 0;
            Int32 i, iNext;
            Boolean found = true;
            Boolean isNegative = (this.xData[t1] < 0.0);
            t2--;
            i = t1;
            while (i < t2)
            {
                iNext = i + 1;
                if (this.xData[i] < 0.0) // assuming that I only have a peak for y > 0
                {// if I cross the x-axis going up
                    isNegative = true;
                }
                if ((this.xData[i] < this.xData[iNext]) && (found))
                {// if I'm climbing up a peak and already found one
                    found = false; // then we should start expecting for the next peak
                }
                else if (this.xData[i] > 0.0)
                {// if I'm at a peak, I only count it if it is greater than 0
                    if ((this.xData[i] > this.xData[iNext]) && (!found) && (isNegative))
                    {// I'm only at a peak if the next point is smaller than the current one, if I haven't already found a peak and if I have previously crossed the x-axis going up
                        found = true;
                        isNegative = false;
                        peakCounter++;
                    }
                }
                i = iNext;
            }
            return peakCounter;
        }
    }
}

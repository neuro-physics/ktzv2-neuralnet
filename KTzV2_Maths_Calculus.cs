using System;

namespace KTzV2.Maths.Calculus
{
    public enum DerivativeMethod
    {
        CentralDiff,
        ImprovedCentralDiff,
        RichardsonExtrap
    }
    public class FirstOrderDerivative
    {
        public Double[] y { get; private set; }
        public Double h { get; private set; }
        private Func<Int32, Double> GetValue_Internal { get; set; }
        public Func<Double> GetError { get; private set; }

        public FirstOrderDerivative(Double[] x, Double[] y, DerivativeMethod m)
            : this(Math.Abs(x[1]-x[0]), y, m)
        {
        }

        public FirstOrderDerivative(Double dx, Double[] y, DerivativeMethod m)
        {
            this.y = y;
            this.h = dx;
            if (m == DerivativeMethod.CentralDiff)
            {
                if (this.y.Length < 3)
                {
                    Console.WriteLine("FirstOrderDerivative: y.Length < 3");
                    //throw new ArgumentException("FirstOrderDerivative: y.Length < 3");
                    this.GetValue_Internal = (x) => Double.NaN;
                    this.GetError = () => Double.NaN;
                }
                else
                {
                    this.GetValue_Internal = this.Central;
                    this.GetError = this.CentralError;
                }
            }
            else if (m == DerivativeMethod.ImprovedCentralDiff)
            {
                if (this.y.Length < 5)
                {
                    Console.WriteLine("FirstOrderDerivative: y.Length < 5");
                    //throw new ArgumentException("FirstOrderDerivative: y.Length < 5");
                    this.GetValue_Internal = (x) => Double.NaN;
                    this.GetError = () => Double.NaN;
                }
                else
                {
                    this.GetValue_Internal = this.Oh4Central;
                    this.GetError = this.Oh4CentralError;
                }
            }
            else if (m == DerivativeMethod.RichardsonExtrap)
            {
                if (this.y.Length < 3)
                {
                    Console.WriteLine("FirstOrderDerivative: y.Length < 3");
                    //throw new ArgumentException("FirstOrderDerivative: y.Length < 3");
                    this.GetValue_Internal = (x) => Double.NaN;
                    this.GetError = () => Double.NaN;
                }
                else
                {
                    this.GetValue_Internal = this.Richardson;
                    this.GetError = this.Oh4CentralError;
                }
            }
        }

        public Double[] GetDerivative()
        {
            Double[] dy = new Double[this.y.Length];
            Int32 i = 0;
            while (i < this.y.Length)
            {
                dy[i] = this.GetValue(i);
                i++;
            }
            return dy;
        }

        public Double[] GetDerivative(out Double[] error)
        {
            Double[] dy = new Double[this.y.Length];
            error = new Double[this.y.Length];
            Int32 i = 0;
            while (i < this.y.Length)
            {
                dy[i] = this.GetValue(i);
                error[i] = this.GetError();
                i++;
            }
            return dy;
        }

        public Double GetValue(Int32 i)
        {
            if ((i < 0) || (i >= this.y.Length))
                return Double.NaN;
            return this.GetValue_Internal(i);
        }

        private Double Central_hHalf(Int32 i)
        {
            if (i < 1)
                return Forward_hHalf(i);
            else if (i > (y.Length - 2))
                return Backward_hHalf(i);
            return (y[i + 1] - y[i - 1]) / h;
        }
        
        private Double Central(Int32 i)
        {
            //if (y == null || y.Length < 3 || i < 1 || i > y.Length - 2 || h == 0)
            if (i < 1)
                return Forward(i);
            else if (i > (y.Length - 2))
                return Backward(i);
            return (y[i + 1] - y[i - 1]) / (2.0D * h);
        }

        private Double Oh4Central(Int32 i)
        {
            //if (y == null || y.Length < 5 || i < 2 || i > y.Length - 3 || h == 0)
            if ((i < 2) || (i > (y.Length - 3)))
                return Central(i);
            return (y[i - 2] - 8.0D * y[i - 1] + 8.0D * y[i + 1] - y[i + 2]) / (12.0D * h);
        }

        private Double Richardson(Int32 i)
        {
            return (4.0D * Central_hHalf(i) - Central(i)) / 3.0D;
        }

        private Double Forward(Int32 i)
        {
            return (-3.0D * y[i] + 4.0D * y[i + 1] - y[i + 2]) / (2.0D * h);
        }

        private Double Backward(Int32 i)
        {
            return (3.0D * y[i] - 4.0D * y[i - 1] + y[i - 2]) / (2.0D * h);
        }

        private Double Forward_hHalf(Int32 i)
        {
            return (-3.0D * y[i] + 4.0D * y[i + 1] - y[i + 2]) / h;
        }

        private Double Backward_hHalf(Int32 i)
        {
            return (3.0D * y[i] - 4.0D * y[i - 1] + y[i - 2]) / h;
        }

        private Double CentralError()
        {
            return this.h * this.h;
        }

        private Double Oh4CentralError()
        {
            return this.h * this.h * this.h * this.h;
        }
    }
}
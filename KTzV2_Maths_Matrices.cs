using System;
using System.Collections.Generic;

namespace KTzV2.Maths.Matrices
{
    /// <summary>
    /// Creates a sparse matrix with elements of the type T (Double, Int, Bool, etc)
    /// </summary>
    /// <typeparam name="T">the type of the elements of the matrix</typeparam>
    public class SparseMatrix<T> where T : struct
    {
        /// <summary>
        /// total number of rows of this matrix
        /// </summary>
        public Int32 nRow { get; private set; }
        /// <summary>
        /// total number of columns of this matrix
        /// </summary>
        public Int32 nCol { get; private set; }

        /// <summary>
        /// a list with all the non-zero elements of this matrix
        /// </summary>
        private List<SparseMatrixElem<T>> elem { get; set; }

        /// <summary>
        /// the index of the next element to get
        /// </summary>
        private Int32 next { get; set; }

        /// <summary>
        /// Initiate a new sparse matrix with the specified number of rows and columns
        /// </summary>
        /// <param name="nRow">number of rows</param>
        /// <param name="nCol">number of columns</param>
        public SparseMatrix(Int32 nRow, Int32 nCol)
        {
            this.nRow = nRow;
            this.nCol = nCol;

            this.elem = new List<SparseMatrixElem<T>>();

            this.ResetNext();
        }

        /// <summary>
        /// Initiate a new sparse matrix with the specified number of rows and columns and the specified elements
        /// </summary>
        /// <param name="nRow">number of rows</param>
        /// <param name="nCol">number of columns</param>
        /// <param name="elems">elements to add</param>
        public SparseMatrix(Int32 nRow, Int32 nCol, List<SparseMatrixElem<T>> elems)
            : this(nRow, nCol)
        {
            this.elem.AddRange(elems);
        }

        /// <summary>
        /// adds an element to the matrix, if it doesn't exists... otherwise, changes the value of the specified element
        /// </summary>
        /// <param name="el">element to be added</param>
        public void AddElement(SparseMatrixElem<T> el)
        {
            // if the element is different from 0 or the default 0 of type T
            // we can add the element... we just ignore the element on the other case
            if (!EqualityComparer<T>.Default.Equals(el.value, default(T)))
            {
                Int32 ind = this.GetElemIndex(el); // checking if the element exists
                if (ind == -1)
                {
                    this.elem.Add(el);
                }
                else
                {
                    this.elem[ind] = el;
                }
            }
        }

        /// <summary>
        /// adds an element without checking if it exists (unsecure)
        /// </summary>
        /// <param name="el">the element</param>
        public void AddElementWithoutCheckingExistence(SparseMatrixElem<T> el)
        {
            // if the element is different from 0 or the default 0 of type T
            // we can add the element... we just ignore the element on the other case
            if (!EqualityComparer<T>.Default.Equals(el.value, default(T)))
            {
                this.elem.Add(el);
            }
        }

        /// <summary>
        /// adds an element with value at the position (i,j)
        /// </summary>
        /// <param name="i">row #</param>
        /// <param name="j">col #</param>
        /// <param name="value">element value</param>
        public void AddElement(Int32 i, Int32 j, T value)
        {
            this.AddElement(new SparseMatrixElem<T>(i, j, value));
        }

        /// <summary>
        /// adds an element without checking if it exists (unsecure)
        /// </summary>
        /// <param name="i">row #</param>
        /// <param name="j">col #</param>
        /// <param name="value">element value</param>
        public void AddElementWithoutCheckingExistence(Int32 i, Int32 j, T value)
        {
            this.AddElementWithoutCheckingExistence(new SparseMatrixElem<T>(i, j, value));
        }

        /// <summary>
        /// removes the element at row i and col j (make it zero, actually)
        /// </summary>
        /// <param name="i">the row # of the element to check</param>
        /// <param name="j">the col # of the element to check</param>
        public void RemElement(Int32 i, Int32 j)
        {
            Int32 ind = this.GetElemIndex(i, j);
            if (ind > 0)
            {
                this.elem.RemoveAt(ind);
                this.elem.TrimExcess();
            }
            else
            {
                throw new Exception("Element doesn't exist! I.e. it is zero.");
            }
        }

        /// <summary>
        /// gets the element at the position (i,j)
        /// </summary>
        /// <param name="i">row #</param>
        /// <param name="j">col #</param>
        /// <returns>element</returns>
        public SparseMatrixElem<T> GetElementAt(Int32 i, Int32 j)
        {
            Int32 ind = this.GetElemIndex(i, j);
            if (ind < 0)
            {
                throw new Exception("Element doesn't exist! I.e. it is zero.");
            }
            return this.elem[ind];
        }

        /// <summary>
        /// gets the element value at the position (i,j)
        /// </summary>
        /// <param name="i">row #</param>
        /// <param name="j">col #</param>
        /// <returns>element value</returns>
        public T GetElementValueAt(Int32 i, Int32 j)
        {
            Int32 ind = this.GetElemIndex(i, j);
            // if the element "doesn't exist", we return the value 0
            // (T)(object)0 is a nasty way of converting to a generic type, but should work for me
            return (ind == -1 ? default(T) : this.elem[ind].value);
        }

        /// <summary>
        /// gets the next non-zero element of the elem list
        /// </summary>
        /// <returns>next non-zero element of the sparse matrix</returns>
        public SparseMatrixElem<T> GetNextElem()
        {
            if (this.next == this.elem.Count - 1)
            {
                this.ResetNext();
            }
            return elem[++this.next];
        }

        /// <summary>
        /// gets the total quantity of non-zero elements on the matrix
        /// </summary>
        /// <returns></returns>
        public Int32 GetElemQuantity()
        {
            return elem.Count;
        }

        /// <summary>
        /// resets the index of the next non-zero element to be returned
        /// </summary>
        private void ResetNext()
        {
            this.next = -1;
        }

        /// <summary>
        /// resets the next parameter to the value val specified
        /// </summary>
        /// <param name="val">the value to reset next to</param>
        private void ResetNext(Int32 val)
        {
            this.next = val;
        }

        /// <summary>
        /// gets the complete matrix represented by this object
        /// </summary>
        /// <returns>a 2d-jagged-array with the matrix</returns>
        public T[][] GetCompleteMatrixByArray()
        {
            Int32 i = 0, j, nextTemp = this.next;
            SparseMatrixElem<T> el;
            this.ResetNext();

            // filling the matrix with 0, assuming that we have most of the elements equal to zero on a sparse matrix
            T[][] matrix = new T[nRow][];
            while (i < this.nRow)
            {
                matrix[i] = new T[nCol];
                j = 0;
                while (j < this.nCol)
                {
                    matrix[i][j] = default(T);
                    j++;
                }
                i++;
            }

            // changing the values of the desired elements
            i = 0;
            j = elem.Count;
            while (i < j)
            {
                el = this.GetNextElem();
                matrix[el.i][el.j] = el.value;
                i++;
            }

            this.ResetNext(nextTemp);
            return matrix;
        }

        /// <summary>
        /// given a column index, j, gets every row index which has non empty entry
        /// </summary>
        /// <param name="j">col #</param>
        /// <param name="elValues">(out) value of element at each row position</param>
        /// <returns>an array with the indeces of the non-empty rows</returns>
        public Int32[] GetNonEmptyRowsInd(Int32 j, out T[] elValues)
        {
            List<SparseMatrixElem<T>> el = elem.FindAll(e => e.j == j);
            Int32 k = 0, n = el.Count;
            elValues = new T[n];
            Int32[] r = new Int32[n];
            while (k < n)
            {
                r[k] = el[k].i;
                elValues[k] = el[k].value;
                k++;
            }
            return r;
        }

        /// <summary>
        /// given a row index, j, gets every col index which has non-empty entry
        /// </summary>
        /// <param name="i">row #</param>
        /// /// <param name="elValues">(out) value of element at each col position</param>
        /// <returns>an array with the indeces of the non-empty cols</returns>
        public Int32[] GetNonEmptyColsInd(Int32 i, out T[] elValues)
        {
            List<SparseMatrixElem<T>> el = elem.FindAll(e => e.i == i);
            Int32 k = 0, n = el.Count;
            elValues = new T[n];
            Int32[] r = new Int32[n];
            while (k < n)
            {
                r[k] = el[k].j;
                elValues[k] = el[k].value;
                k++;
            }
            return r;
        }

        /// <summary>
        /// verifies if the specified element exists... doesn't matter the value of the element
        /// </summary>
        /// <param name="el">element to search for</param>
        /// <returns>the index of the specified element, if it's found... otherwise it returns -1</returns>
        private Int32 GetElemIndex(SparseMatrixElem<T> el)
        {
            return this.elem.FindIndex(e => e.i == el.i && e.j == el.j);
        }

        /// <summary>
        /// verifies if the specified element exists... doesn't matter the value of the element
        /// </summary>
        /// <param name="i">the row # of the element to check</param>
        /// <param name="j">the col # of the element to check</param>
        /// <returns>the index of the specified element, if it's found... otherwise it returns -1</returns>
        private Int32 GetElemIndex(Int32 i, Int32 j)
        {
            return this.elem.FindIndex(e => e.i == i && e.j == j);
        }
    }

    /// <summary>
    /// a sparse matrix element
    /// </summary>
    /// <typeparam name="T">the type of the element</typeparam>
    public class SparseMatrixElem<T> where T : struct
    {
        public Int32 i { get; set; }
        public Int32 j { get; set; }
        public T value { get; set; }

        /// <summary>
        /// creates an sparse matrix element at the specifieds row and col, with determined value of type T
        /// </summary>
        /// <param name="i">the row number of the element</param>
        /// <param name="j">the col number of the element</param>
        /// <param name="value">the value of the element</param>
        public SparseMatrixElem(Int32 i, Int32 j, T value)
        {
            this.i = i;
            this.j = j;
            this.value = value;
        }
    }
}
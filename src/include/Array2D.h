#ifndef ARRAY2D_H
#define ARRAY2D_H

#include <iostream>
#include <vector>

template <class T>
class Array2D
{
private:
    std::vector<T> StridedSlice( int start, int length, int stride ) const;

    int m_rows;
    int m_cols;

    std::vector<T> m_data;

public:
    Array2D();
    Array2D( int rows, int cols, const T& initVal = T() );
    Array2D( const Array2D<T>& );

    // Size and structure
    void Resize( int rows, int cols );
    int NumRows() const                       { return m_rows; }
    int NumCols() const                    { return m_cols; }
    int Size() const                   { return m_data.size(); }

    // Direct vector access and indexing
    //operator const std::vector<T>& () const        { return m_data; }
    const std::vector<T>& DataBlock() const	{ return m_data; }
    int Index( int row, int col ) const		{ return row * m_cols + col; }

    // Get a single value
    //      T & Value( int row, int col )       { return m_data[Index(row,col)]; }
    //const T & Value( int row, int col ) const { return m_data[Index(row,col)]; }
          T & operator[]( size_t idx )        { return m_data[idx]; }
    const T & operator[]( size_t idx ) const  { return m_data[idx]; }
	  T& operator() (const int row, const int col)		{ return m_data[Index(row,col)]; }
    const T& operator() (const int row, const int col) const	{ return m_data[Index(row,col)]; }

    /* copy the contents of another Array2D into this one as a sub matrix */
    bool CopySub(const Array2D&, int rowbeg, int colbeg ); 
    // Simple row or column slices
    std::vector<T> GetRow( int row, int colBegin = 0, int colEnd = -1 ) const;
    std::vector<T> GetCol( int row, int colBegin = 0, int colEnd = -1 ) const;

};


template <class T>
std::vector<T> Array2D<T>::StridedSlice( int start, int length, int stride ) const
{
    std::vector<T> result;
    result.reserve( length );
    const T *pos = &m_data[start];
    for( int i = 0; i < length; i++ ) {
        result.push_back(*pos);
        pos += stride;
    }
    return result;
}

template <class T>
Array2D<T>::Array2D() {}

template <class T>
Array2D<T>::Array2D( int rows, int cols, const T& initVal )
    : m_data( rows * cols, initVal )
    , m_rows( rows )
    , m_cols( cols ) {}


template <class T>
Array2D<T>::Array2D( const Array2D<T>& Ain ) 
   : m_data( Ain.DataBlock() )
   , m_rows( Ain.NumRows() )
   , m_cols( Ain.NumCols() )
{}


template <class T>
void Array2D<T>::Resize(int rows, int cols ) {
   m_rows = rows;
   m_cols = cols;
   m_data.resize( rows * cols );
}


template <class T>
bool Array2D<T>::CopySub(const Array2D& Ain, int rowbeg, int colbeg ) {
   if( m_rows-rowbeg < Ain.NumRows() || m_cols-colbeg < Ain.NumCols() ) return false;
   int irow = rowbeg, icol;
   for( int irowin=0; irowin<Ain.NumRows(); irow++, irowin++ ) {
      icol = colbeg;
      for( int icolin=0; icolin<Ain.NumCols(); icol++, icolin++ ) (*this)(irow, icol) = Ain(irowin, icolin);
   }
   return true;
}

/* construct an vector storing a single row or column with the StridedSlice function */
template <class T>
std::vector<T> Array2D<T>::GetRow( int row, int colBegin, int colEnd ) const {
    if( colEnd < 0 ) colEnd = m_cols-1;
    if( colBegin <= colEnd )
        return StridedSlice( Index(row,colBegin), colEnd-colBegin+1, 1 );
    else
        return StridedSlice( Index(row,colBegin), colBegin-colEnd+1, -1 );
}

template <class T>
std::vector<T> Array2D<T>::GetCol( int col, int rowBegin, int rowEnd ) const {
    if( rowEnd < 0 ) rowEnd = m_rows-1;
    if( rowBegin <= rowEnd )
        return StridedSlice( Index(rowBegin,col), rowEnd-rowBegin+1, m_cols );
    else
        return StridedSlice( Index(rowBegin,col), rowBegin-rowEnd+1, -m_cols );
}

#endif

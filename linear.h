/*
 linear
 
 Copyright (c) 2010 Cass Everitt
 
 from glh -
 
 Copyright (c) 2000-2009 Cass Everitt
 Copyright (c) 2000 NVIDIA Corporation
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or
 without modification, are permitted provided that the following
 conditions are met:
 
 * Redistributions of source code must retain the above
 copyright notice, this list of conditions and the following
 disclaimer.
 
 * Redistributions in binary form must reproduce the above
 copyright notice, this list of conditions and the following
 disclaimer in the documentation and/or other materials
 provided with the distribution.
 
 * The names of contributors to this software may not be used
 to endorse or promote products derived from this software
 without specific prior written permission. 
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 REGENTS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 POSSIBILITY OF SUCH DAMAGE. 
 
 
 Cass Everitt
 */

// Author:  Cass W. Everitt

#ifndef __R3_LINEAR_H__
#define __R3_LINEAR_H__

//#include <memory.h>
#include <math.h>
#include <assert.h>
#include <algorithm>

#ifdef _WIN32
# define TEMPLATE_FUNCTION
#else
# define TEMPLATE_FUNCTION <>
#endif

#define     R3_QUATERNION_NORMALIZATION_THRESHOLD  64

#define     R3_RAD_TO_DEG      float(57.2957795130823208767981548141052)
#define     R3_DEG_TO_RAD      float(0.0174532925199432957692369076848861)
#define     R3_ZERO            float(0.0)
#define     R3_ONE             float(1.0)
#define     R3_TWO             float(2.0)
#define     R3_EPSILON         float(10e-8)
#define     R3_PI              float(3.1415926535897932384626433832795)    

#define     r3_equivalent(a,b)     (((a < b + R3_EPSILON) && (a > b - R3_EPSILON)) ? true : false)
#define     r3_greater_than(a,b)   ( (a > b - R3_EPSILON) ? true : false)
#define     r3_less_than(a,b)      ( (a < b + R3_EPSILON) ? true : false)

namespace r3 {
	
	inline float ToDegrees(float radians) { 
		return radians * R3_RAD_TO_DEG; 
	}
	
	inline float ToRadians(float degrees) { 
		return degrees * R3_DEG_TO_RAD;
	}
	
	template < typename T > class Line;
	template < typename T > class Plane;
	template < typename T > class Matrix3;
	template < typename T > class Matrix4;
	template < typename T > class Quaternion;

	template < class T >
	struct Vec2  {

		Vec2() { x = y = 0.f; }
		Vec2(const T * tp) { x = tp[0]; y = tp[1]; }

		Vec2(T x, T y ) {
			Ptr()[0] = x;
			Ptr()[1] = y; 
		}
		
		void GetValue(T & x, T & y) const {
			x = Ptr()[0]; y = Ptr()[1]; 
		}
		
		Vec2 & SetValue( const T & x, const T & y) {
			Ptr()[0] = x; Ptr()[1] = y; return *this;
		}
		
		T x, y;

		T * Ptr() { return &x; }
		const T * Ptr() const { return &x; }
		
		static const int N = 2;
		// generic part
		
		int Size() const {
			return N;
		}
				
		const T * GetValue() const {
			return Ptr(); 
		}
		
		
		T Dot( const Vec2 & rhs ) const { 
			T r = 0;
			for(int i = 0; i < N; i++) r += Ptr()[i]*rhs.Ptr()[i];
			return r;
		}
		
		T Length() const {
			T r = 0;
			for(int i = 0; i < N; i++) r += Ptr()[i]*Ptr()[i]; 
			return T(sqrt(r));
		}	
		
		T SquareNorm() const {
			T r = 0;
			for(int i = 0; i < N; i++) r += Ptr()[i]*Ptr()[i]; 
			return r;
		}	
		
		void  Negate() {
			for(int i = 0; i < N; i++) Ptr()[i] = -Ptr()[i]; 
		}
		
		
		T Normalize() { 
			T sum(0);
			for(int i = 0; i < N; i++) 
				sum += Ptr()[i]*Ptr()[i];
			sum = T(sqrt(sum));
			if (sum > R3_EPSILON)
				for(int i = 0; i < N; i++) 
					Ptr()[i] /= sum;
			return sum;
		}
		
		
		Vec2 & SetValue( const T * rhs ) {
			for(int i = 0; i < N; i++) Ptr()[i] = rhs[i]; return *this; 
		}
		
		T & operator [] ( int i ) {
			return Ptr()[i]; 
		}
		
		const T & operator [] ( int i ) const {
			return Ptr()[i]; 
		}
		
		Vec2 & operator *= ( T d ) {
			for(int i = 0; i < N; i++) Ptr()[i] *= d; return *this;
		}
		
		Vec2 & operator *= ( const Vec2 & u ) {
			for(int i = 0; i < N; i++) Ptr()[i] *= u[i]; return *this;
		}
		
		Vec2 & operator /= ( T d ) {
			if(d == 0) return *this; for(int i = 0; i < N; i++) Ptr()[i] /= d; return *this;
		}
		
		Vec2 & operator += ( const Vec2 & u ) {
			for(int i = 0; i < N; i++) Ptr()[i] += u.Ptr()[i]; return *this;
		}
		
		Vec2 & operator -= ( const Vec2 & u ) {
			for(int i = 0; i < N; i++) Ptr()[i] -= u.Ptr()[i]; return *this;
		}
		
		
		Vec2 operator - () const {
			Vec2 rv(*this); rv.Negate(); return rv; 
		}
		
		Vec2 operator + ( const Vec2 &rhs) const {
			Vec2 rt(*this); return rt += rhs; 
		}
		
		Vec2 operator - ( const Vec2 &rhs) const {
			Vec2 rt(*this); return rt -= rhs; 
		}
		
	};

	
	// vector friend operators
	
	template <class T> inline
	Vec2<T> operator * ( const Vec2<T> & b, T d ) {
		Vec2<T> rt(b);
		return rt *= d;
	}
	
	template <class T> inline
	Vec2<T> operator * ( T d, const Vec2<T> & b ) {
		return b*d; 
	}
	
	template <class T> inline
	Vec2<T> operator * ( const Vec2<T> & b, const Vec2<T> & d )	{
		Vec2<T> rt(b);
		return rt *= d;
	}
	
	template <class T> inline
	Vec2<T> operator / ( const Vec2<T> & b, T d ) {
		Vec2<T> rt(b); return rt /= d; 
	}
	
	template <class T> inline 
	Vec2<T> operator + ( const Vec2<T> & v1, const Vec2<T> & v2 ) {
		Vec2<T> rt(v1); return rt += v2; 
	}
	
	template <class T> inline
	Vec2<T> operator - ( const Vec2<T> & v1, const Vec2<T> & v2 ) {
		Vec2<T> rt(v1); return rt -= v2; 
	}
	
	
	template <class T> inline
	bool operator == ( const Vec2<T> & v1, const Vec2<T> & v2 ) {
		for(int i = 0; i < 2; i++)
			if(v1.Ptr()[i] != v2.Ptr()[i])
				return false;
		return true;
	}
	
	template <class T> inline
	bool operator != ( const Vec2<T> & v1, const Vec2<T> & v2 ) {
		return !(v1 == v2); 
	}

	template <typename T> inline
		Vec2<T> Min( const Vec2<T> & v1, const Vec2<T> &v2 ) {
			Vec2<T> r;
			r.x = std::min( v1.x, v2.x );
			r.y = std::min( v1.y, v2.y );
			return r;
	}

	template <typename T> inline
		Vec2<T> Max( const Vec2<T> & v1, const Vec2<T> &v2 ) {
			Vec2<T> r;
			r.x = std::max( v1.x, v2.x );
			r.y = std::max( v1.y, v2.y );
			return r;
	}
	
	
	
	template < typename T>
	class Vec3 {
    public:
		Vec3() { x = y = z = 0.f; }
		Vec3(const T * tp) { x = tp[0]; y = tp[1]; z = tp[2]; }
		
		Vec3(T x, T y, T z) {
			Ptr()[0] = x; Ptr()[1] = y; Ptr()[2] = z; 
		}
		
		void GetValue(T & x, T & y, T & z) const {
			x = Ptr()[0]; y = Ptr()[1]; z = Ptr()[2]; 
		}
		
		Vec3 Cross( const Vec3 &rhs ) const	{
			Vec3 rt;
			rt.x = Ptr()[1] * rhs.z - Ptr()[2] * rhs.y;
			rt.y = Ptr()[2] * rhs.x - Ptr()[0] * rhs.z;
			rt.z = Ptr()[0] * rhs.y - Ptr()[1] * rhs.x;	
			return rt;
		}
		
		Vec3 & SetValue( const T & x, const T & y, const T & z) {
			Ptr()[0] = x; Ptr()[1] = y; Ptr()[2] = z; return *this; 
		}
		
		T x, y, z;
		
		T * Ptr() { return &x; }
		const T * Ptr() const { return &x; }
		
		static const int N = 3;
		// generic part
		
		int Size() const {
			return N;
		}
				
		const T * GetValue() const {
			return Ptr(); 
		}
		
		
		T Dot( const Vec3 & rhs ) const { 
			T r = 0;
			for(int i = 0; i < N; i++) r += Ptr()[i]*rhs.Ptr()[i];
			return r;
		}
		
		T Length() const {
			T r = 0;
			for(int i = 0; i < N; i++) r += Ptr()[i]*Ptr()[i]; 
			return T(sqrt(r));
		}	
		
		T SquareNorm() const {
			T r = 0;
			for(int i = 0; i < N; i++) r += Ptr()[i]*Ptr()[i]; 
			return r;
		}	
		
		void  Negate() {
			for(int i = 0; i < N; i++) Ptr()[i] = -Ptr()[i]; 
		}
		
		
		T Normalize() { 
			T sum(0);
			for(int i = 0; i < N; i++) 
				sum += Ptr()[i]*Ptr()[i];
			sum = T(sqrt(sum));
			if (sum > R3_EPSILON)
				for(int i = 0; i < N; i++) 
					Ptr()[i] /= sum;
			return sum;
		}
		
		
		Vec3 & SetValue( const T * rhs ) {
			for(int i = 0; i < N; i++) Ptr()[i] = rhs[i]; return *this; 
		}
		
		T & operator [] ( int i ) {
			return Ptr()[i]; 
		}
		
		const T & operator [] ( int i ) const {
			return Ptr()[i]; 
		}
		
		Vec3 & operator *= ( T d ) {
			for(int i = 0; i < N; i++) Ptr()[i] *= d; return *this;
		}
		
		Vec3 & operator *= ( const Vec3 & u ) {
			for(int i = 0; i < N; i++) Ptr()[i] *= u[i]; return *this;
		}
		
		Vec3 & operator /= ( T d ) {
			if(d == 0) return *this; for(int i = 0; i < N; i++) Ptr()[i] /= d; return *this;
		}
		
		Vec3 & operator += ( T d ) {
			for(int i = 0; i < N; i++) Ptr()[i] += d; return *this;
		}

		Vec3 & operator += ( const Vec3 & u ) {
			for(int i = 0; i < N; i++) Ptr()[i] += u.Ptr()[i]; return *this;
		}
		
		Vec3 & operator -= ( T d ) {
			for(int i = 0; i < N; i++) Ptr()[i] -= d; return *this;
		}
		
		Vec3 & operator -= ( const Vec3 & u ) {
			for(int i = 0; i < N; i++) Ptr()[i] -= u.Ptr()[i]; return *this;
		}
		
		
		Vec3 operator - () const {
			Vec3 rv(*this); rv.Negate(); return rv; 
		}
		
		Vec3 operator + ( const Vec3 &rhs) const {
			Vec3 rt(*this); return rt += rhs; 
		}
		
		Vec3 operator - ( const Vec3 &rhs) const {
			Vec3 rt(*this); return rt -= rhs; 
		}
		
	};
	
	
	// vector friend operators
	
	template <class T> inline
	Vec3<T> operator * ( const Vec3<T> & b, T d ) {
		Vec3<T> rt(b);
		return rt *= d;
	}
	
	template <class T> inline
	Vec3<T> operator * ( T d, const Vec3<T> & b ) {
		return b*d; 
	}
	
	template <class T> inline
	Vec3<T> operator * ( const Vec3<T> & b, const Vec3<T> & d )	{
		Vec3<T> rt(b);
		return rt *= d;
	}
	
	template <class T> inline
	Vec3<T> operator / ( const Vec3<T> & b, T d ) {
		Vec3<T> rt(b); return rt /= d; 
	}
	
	template <class T> inline 
	Vec3<T> operator + ( const Vec3<T> & v1, const Vec3<T> & v2 ) {
		Vec3<T> rt(v1); return rt += v2; 
	}
	
	template <class T> inline
	Vec3<T> operator - ( const Vec3<T> & v1, const Vec3<T> & v2 ) {
		Vec3<T> rt(v1); return rt -= v2; 
	}
	
	
	template <class T> inline
	bool operator == ( const Vec3<T> & v1, const Vec3<T> & v2 ) {
		for(int i = 0; i < 3; i++)
			if(v1.Ptr()[i] != v2.Ptr()[i])
				return false;
		return true;
	}
	
	template <class T> inline
	bool operator != ( const Vec3<T> & v1, const Vec3<T> & v2 ) {
		return !(v1 == v2); 
	}
	
	template <typename T> inline
		Vec3<T> Min( const Vec3<T> & v1, const Vec3<T> &v2 ) {
			Vec3<T> r;
			r.x = std::min( v1.x, v2.x );
			r.y = std::min( v1.y, v2.y );
			r.z = std::min( v1.z, v2.z );
			return r;
	}

	template <typename T> inline
		Vec3<T> Max( const Vec3<T> & v1, const Vec3<T> &v2 ) {
			Vec3<T> r;
			r.x = std::max( v1.x, v2.x );
			r.y = std::max( v1.y, v2.y );
			r.z = std::max( v1.z, v2.z );
			return r;
	}

	
	
	
	
	
	
	template < typename T>
    class Vec4 {
    public:
		Vec4() { x = y = z = 0.f; w = 1.f; }
		Vec4(const T * tp) { x = tp[0]; y = tp[1]; z = tp[2]; w = tp[3]; }
        Vec4(const Vec3<T> & t, T fourth) {
			Ptr()[0] = t.x; Ptr()[1] = t.y; Ptr()[2] = t.z; Ptr()[3] = fourth; 
		}
        Vec4(T x, T y, T z = 0, T w = 1) {
			Ptr()[0] = x; Ptr()[1] = y; Ptr()[2] = z; Ptr()[3] = w;
		}
		
        void GetValue(T & x, T & y, T & z, T & w) const {
			x = Ptr()[0]; y = Ptr()[1]; z = Ptr()[2]; w = Ptr()[3];
		}
		
        Vec4 & SetValue( const T & x, const T & y, const T & z, const T & w) {
			Ptr()[0] = x; Ptr()[1] = y; Ptr()[2] = z; Ptr()[3] = w; return *this;
		}

		T x, y, z, w;
		
		T * Ptr() { return &x; }
		const T * Ptr() const { return &x; }
		
		static const int N = 4;
		// generic part
		
		int Size() const {
			return N;
		}
		
		
		const T * GetValue() const {
			return Ptr(); 
		}
		
		
		T Dot( const Vec4 & rhs ) const { 
			T r = 0;
			for(int i = 0; i < N; i++) r += Ptr()[i]*rhs.Ptr()[i];
			return r;
		}
		
		T Length() const {
			T r = 0;
			for(int i = 0; i < N; i++) r += Ptr()[i]*Ptr()[i]; 
			return T(sqrt(r));
		}	
		
		T SquareNorm() const {
			T r = 0;
			for(int i = 0; i < N; i++) r += Ptr()[i]*Ptr()[i]; 
			return r;
		}	
		
		void  Negate() {
			for(int i = 0; i < N; i++) Ptr()[i] = -Ptr()[i]; 
		}
		
		
		T Normalize() { 
			T sum(0);
			for(int i = 0; i < N; i++) 
				sum += Ptr()[i]*Ptr()[i];
			sum = T(sqrt(sum));
			if (sum > R3_EPSILON)
				for(int i = 0; i < N; i++) 
					Ptr()[i] /= sum;
			return sum;
		}
		
		
		Vec4 & SetValue( const T * rhs ) {
			for(int i = 0; i < N; i++) Ptr()[i] = rhs[i]; return *this; 
		}
		
		T & operator [] ( int i ) {
			return Ptr()[i]; 
		}
		
		const T & operator [] ( int i ) const {
			return Ptr()[i]; 
		}
		
		Vec4 & operator *= ( T d ) {
			for(int i = 0; i < N; i++) Ptr()[i] *= d; return *this;
		}
		
		Vec4 & operator *= ( const Vec4 & u ) {
			for(int i = 0; i < N; i++) Ptr()[i] *= u[i]; return *this;
		}
		
		Vec4 & operator /= ( T d ) {
			if(d == 0) return *this; for(int i = 0; i < N; i++) Ptr()[i] /= d; return *this;
		}
		
		Vec4 & operator += ( const Vec4 & u ) {
			for(int i = 0; i < N; i++) Ptr()[i] += u.Ptr()[i]; return *this;
		}
		
		Vec4 & operator -= ( const Vec4 & u ) {
			for(int i = 0; i < N; i++) Ptr()[i] -= u.Ptr()[i]; return *this;
		}
		
		
		Vec4 operator - () const {
			Vec4 rv(*this); rv.negate(); return rv; 
		}
		
		Vec4 operator + ( const Vec4 &rhs) const {
			Vec4 rt(*this); return rt += rhs; 
		}
		
		Vec4 operator - ( const Vec4 &rhs) const {
			Vec4 rt(*this); return rt -= rhs; 
		}
		

	};
	
    template < typename T > inline
    Vec3<T> Homogenize(const Vec4<T> & v) {
		Vec3<T> rt;
		assert(v[3] != R3_ZERO);
		rt[0] = v[0]/v[3];
		rt[1] = v[1]/v[3];
		rt[2] = v[2]/v[3];
		return rt;
    }
	
	// vector friend operators
	
	template <class T> inline
	Vec4<T> operator * ( const Vec4<T> & b, T d ) {
		Vec4<T> rt(b);
		return rt *= d;
	}
	
	template <class T> inline
	Vec4<T> operator * ( T d, const Vec4<T> & b ) {
		return b*d; 
	}
	
	template <class T> inline
	Vec4<T> operator * ( const Vec4<T> & b, const Vec4<T> & d )	{
		Vec4<T> rt(b);
		return rt *= d;
	}
	
	template <class T> inline
	Vec4<T> operator / ( const Vec4<T> & b, T d ) {
		Vec4<T> rt(b); return rt /= d; 
	}
	
	template <class T> inline 
	Vec4<T> operator + ( const Vec4<T> & v1, const Vec4<T> & v2 ) {
		Vec4<T> rt(v1); return rt += v2; 
	}
	
	template <class T> inline
	Vec4<T> operator - ( const Vec4<T> & v1, const Vec4<T> & v2 ) {
		Vec4<T> rt(v1); return rt -= v2; 
	}
	
	
	template <class T> inline
	bool operator == ( const Vec4<T> & v1, const Vec4<T> & v2 ) {
		for(int i = 0; i < 4; i++)
			if(v1.Ptr()[i] != v2.Ptr()[i])
				return false;
		return true;
	}
	
	template <class T> inline
	bool operator != ( const Vec4<T> & v1, const Vec4<T> & v2 ) {
		return !(v1 == v2); 
	}
	
	template <typename T> inline
		Vec4<T> Min( const Vec4<T> & v1, const Vec4<T> &v2 ) {
			Vec4<T> r;
			r.x = std::min( v1.x, v2.x );
			r.y = std::min( v1.y, v2.y );
			r.z = std::min( v1.z, v2.z );
			r.w = std::min( v1.w, v2.w );
			return r;
	}

	template <typename T> inline
		Vec4<T> Max( const Vec4<T> & v1, const Vec4<T> &v2 ) {
			Vec4<T> r;
			r.x = std::max( v1.x, v2.x );
			r.y = std::max( v1.y, v2.y );
			r.z = std::max( v1.z, v2.z );
			r.w = std::max( v1.w, v2.w );
			return r;
	}


	
	
	
	
	
	template < typename T>
    class Line	{
	public:
		
		Line() {
			SetValue(Vec3<T>(0,0,0),Vec3<T>(0,0,1)); 
		}
		
		Line( const Vec3<T> & p0, const Vec3<T> &p1) {
			SetValue(p0,p1); 
		}
		
		void SetValue( const Vec3<T> &p0, const Vec3<T> &p1) {
			position = p0;
			direction = p1-p0;
			direction.normalize();
		}
		
		bool GetClosestPoints(const Line &line2, Vec3<T> &pointOnThis, Vec3<T> &pointOnThat) {
			
			// quick check to see if parallel -- if so, quit.
			if(fabs(direction.dot(line2.direction)) == 1.0)
				return 0;
			Line l2 = line2;
			
			// Algorithm: Brian Jean
			// 
			register T u;
			register T v;
			Vec3<T> Vr = direction;
			Vec3<T> Vs = l2.direction;
			register T Vr_Dot_Vs = Vr.dot(Vs);
			register T detA = T(1.0 - (Vr_Dot_Vs * Vr_Dot_Vs));
			Vec3<T> C = l2.position - position;
			register T C_Dot_Vr =  C.dot(Vr);
			register T C_Dot_Vs =  C.dot(Vs);
			
			u = (C_Dot_Vr - Vr_Dot_Vs * C_Dot_Vs)/detA;
			v = (C_Dot_Vr * Vr_Dot_Vs - C_Dot_Vs)/detA;
			
			pointOnThis = position;
			pointOnThis += direction * u;
			pointOnThat = l2.position;
			pointOnThat += l2.direction * v;
			
			return 1;
		}
		
		Vec3<T> GetClosestPoint(const Vec3<T> &point) {
			Vec3<T> np = point - position;
			Vec3<T> rp = direction*direction.dot(np)+position;
			return rp;
		}
		
		const Vec3<T> & GetPosition() const {return position;}
		
		const Vec3<T> & GetDirection() const {return direction;}
		
		//protected:
		Vec3<T> position;
		Vec3<T> direction;
	};
	
	
	template <typename T>
	struct LineSegment2 {
		Vec2<T> a, b;
		LineSegment2() {}
		LineSegment2( const Vec2<T> & ptA, const Vec2<T> & ptB ) 
		: a( ptA ), b( ptB ) {}
		Vec3<T> GetPlane() const { // not normalized
			Vec3<T> p;
			p.x = a.y - b.y;
			p.y = b.x - a.x;
			p.z = - ( p.x * a.x + p.y * a.y );
			return p;
		}
	};
	
	template <typename T>
	inline bool Intersect( const LineSegment2<T> & s0, const LineSegment2<T> & s1 ) {
		Vec3<T> p = s0.GetPlane();
		if ( p.Dot( Vec3<T>( s1.a.x, s1.a.y, 1 ) ) * p.Dot( Vec3<T>( s1.b.x, s1.b.y, 1 ) ) > 0 ) {
			return false;
		}
		p = s1.GetPlane();
		if ( p.Dot( Vec3<T>( s0.a.x, s0.a.y, 1 ) ) * p.Dot( Vec3<T>( s0.b.x, s0.b.y, 1 ) ) > 0 ) {
			return false;
		}
		return true;
	}
	
	
	
	
	
	
	
	
	
	// Matrix3 
	template <typename T>
	class Matrix3 {
	public:
		T m[3][3];
		
		Matrix3() {
			MakeIdentity();
		}
    
		T & element (int row, int col) {
			return m[row][col];
		}
		
		const T & element (int row, int col) const {
			return m[row][col];
		}

    void SetValue( const T & t ) {
      for( int i = 0; i < 3; i++ ) {
        for( int j = 0; j < 3; j++ ) {
          element( i, j ) = t;
        }
      }
    }
    
		
		template <typename TIN>
		Matrix3( const TIN * in ) {
			m[0][0] = in[0];  m[0][1] = in[1];  m[0][2] = in[2];
			m[1][0] = in[3];  m[1][1] = in[4];  m[1][2] = in[5];
			m[2][0] = in[6];  m[2][1] = in[7];  m[2][2] = in[8];
		}
		
		void MakeIdentity() {
			m[0][0] = 1;  m[0][1] = 0;  m[0][2] = 0;
			m[1][0] = 0;  m[1][1] = 1;  m[1][2] = 0;
			m[2][0] = 0;  m[2][1] = 0;  m[2][2] = 1;
		}
		
		Matrix3 Adjugate() const { 
			Matrix3 m3;
			int L[3] = { 1, 0, 0 };
			int G[3] = { 2, 2, 1 };
			for ( int row = 0; row < 3; row++ ) {
				for ( int col = 0; col < 3; col++ ) {
					T t
					= m[ L[row] ][ L[col] ] * m[ G[row] ][ G[col] ]
					- m[ G[row] ][ L[col] ] * m[ L[row] ][ G[col] ];
					
					m3.m[row][col] = ( (row + col ) & 0x1 ) ? -t : t;
				}
			}
			return m3.Transpose();
		}
		
		Matrix3 Transpose() const {
			Matrix3 m3;
			m3.m[0][0] = m[0][0]; m3.m[1][0] = m[0][1]; m3.m[2][0] = m[0][2];
			m3.m[0][1] = m[1][0]; m3.m[1][1] = m[1][1]; m3.m[2][1] = m[1][2];
			m3.m[0][2] = m[2][0]; m3.m[1][2] = m[2][1]; m3.m[2][2] = m[2][2];
			return m3;
		}

    Matrix3 & MultRight( const Matrix3 & b ) {
			Matrix3 mt(*this);
			SetValue(T(0));
			
			for(int i=0; i < 3; i++)
				for(int j=0; j < 3; j++)
					for(int c=0; c < 3; c++)
						element(i,j) += mt(i,c) * b(c,j);
			return *this;
		}
		
		Matrix3 & MultLeft( const Matrix3 & b ) {
			Matrix3 mt(*this);
			SetValue(T(0));
			
			for(int i=0; i < 3; i++)
				for(int j=0; j < 3; j++)
					for(int c=0; c < 3; c++)
						element(i,j) += b(i,c) * mt(c,j);
			return *this;
		}
		
		T Determinant() const {
			T result
			= m[0][0] * m[1][1] * m[2][2]
			+ m[0][1] * m[1][2] * m[2][0]
			+ m[0][2] * m[1][0] * m[2][1]
			- m[2][0] * m[1][1] * m[0][2]
			- m[2][1] * m[1][2] * m[0][0]
			- m[2][2] * m[1][0] * m[0][1];
			return result;
		}
		
		void Div( T t ) {
			Mul( T(1.0) / t );
		}
		
		void Mul( T t ) {
			for ( int row = 0; row < 3; row++ ) {
				for ( int col = 0; col < 3; col++ ) {
					m[row][col] *= t;
				}
			}
		}
		
		// Cramer's Rule
		Matrix3 Inverse() const {
			Matrix3 m3 = Adjugate();
			m3.Div( Determinant() );
			return m3;
		}
		
		Vec3<T> GetRow( int i ) const {
			return Vec3<T>( m[i][0], m[i][1], m[i][2] );			
		}
		
		void SetRow( int i, const Vec3<T> &v ) {
			m[i][0] = v.x;  m[i][1] = v.y;  m[i][2] = v.z;
		}
		
		Vec3<T> GetColumn( int i ) const {
			return Vec3<T>( m[0][i], m[1][i], m[2][i] );			
		}
		
		void SetColumn( int i, const Vec3<T> &v ) {
			m[0][i] = v.x;  m[1][i] = v.y;  m[2][i] = v.z;
		}
		
		T & operator() ( int row, int col ) {
			return m[row][col];
		}
		
		const T & operator() ( int row, int col ) const {
			return m[row][col];
		}
	};
	
	template <typename T> inline
	Vec3<T> operator * ( const Matrix3<T> &m, const Vec3<T> &v) {
		return Vec3<T>(m(0,0) * v.x + m(0,1) * v.y + m(0,2) * v.z,
					   m(1,0) * v.x + m(1,1) * v.y + m(1,2) * v.z,
					   m(2,0) * v.x + m(2,1) * v.y + m(2,2) * v.z);	
	}
	
  template <typename T> inline
  Matrix3<T> operator * ( const Matrix3<T> & m1, const Matrix3<T> & m2 ) {
    static Matrix3<T> product;
      
    product = m1;
    product.MultRight(m2);
      
    return product;
  }
    

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// Matrix4
	
	template < typename T>
	class Matrix4 {
		
	public:
		
		Matrix4() { MakeIdentity(); }
		
		Matrix4( T * m ) { SetValue(m); }
		
		Matrix4( T a00, T a01, T a02, T a03,
				T a10, T a11, T a12, T a13,
				T a20, T a21, T a22, T a23,
				T a30, T a31, T a32, T a33 ) {
			element(0,0) = a00;
			element(0,1) = a01;
			element(0,2) = a02;
			element(0,3) = a03;
			
			element(1,0) = a10;
			element(1,1) = a11;
			element(1,2) = a12;
			element(1,3) = a13;
			
			element(2,0) = a20;
			element(2,1) = a21;
			element(2,2) = a22;
			element(2,3) = a23;
			
			element(3,0) = a30;
			element(3,1) = a31;
			element(3,2) = a32;
			element(3,3) = a33;
		}
		
		
		void GetValue( T * mp ) const {
			int c = 0;
			for(int j=0; j < 4; j++)
				for(int i=0; i < 4; i++)
					mp[c++] = element(i,j);
		}
		
		
		const T * GetValue() const	{
			return m;
		}
		
		void SetValue( T * mp) {
			int c = 0;
			for(int j=0; j < 4; j++)
				for(int i=0; i < 4; i++)
					element(i,j) = mp[c++];
		}
		
		void SetValue( T r ) {
			for(int i=0; i < 4; i++)
				for(int j=0; j < 4; j++)
					element(i,j) = r;
		}
		
		void MakeIdentity() {
			element(0,0) = 1.0;
			element(0,1) = 0.0;
			element(0,2) = 0.0; 
			element(0,3) = 0.0;
			
			element(1,0) = 0.0;
			element(1,1) = 1.0; 
			element(1,2) = 0.0;
			element(1,3) = 0.0;
			
			element(2,0) = 0.0;
			element(2,1) = 0.0;
			element(2,2) = 1.0;
			element(2,3) = 0.0;
			
			element(3,0) = 0.0; 
			element(3,1) = 0.0; 
			element(3,2) = 0.0;
			element(3,3) = 1.0;
		}
		
		
		static Matrix4 Identity() {
			static Matrix4 mident (
								   1.0, 0.0, 0.0, 0.0,
								   0.0, 1.0, 0.0, 0.0,
								   0.0, 0.0, 1.0, 0.0,
								   0.0, 0.0, 0.0, 1.0  );
			return mident;
		}
		
		static Matrix4 Scale( T s ) {
			Matrix4 m;
			m.SetScale( s );
			return m;
		}
		
		static Matrix4 Scale( Vec3<T> s ) {
			Matrix4 m;
			m.SetScale( s );
			return m;
		}

		static Matrix4 Translate( Vec3<T> t ) {
			Matrix4 m;
			m.SetTranslate( t );
			return m;
		}
		
		void SetScale( T s ) {
			element(0,0) = s;
			element(1,1) = s;
			element(2,2) = s;
		}
		
		void SetScale( const Vec3<T> & s ) {
			element(0,0) = s.x;
			element(1,1) = s.y;
			element(2,2) = s.z;
		}
		
		
		void SetTranslate( const Vec3<T> & t )	{
			element(0,3) = t.x;
			element(1,3) = t.y;
			element(2,3) = t.z;
		}
		
		void SetRow(int r, const Vec4<T> & t) {
			element(r,0) = t.x;
			element(r,1) = t.y;
			element(r,2) = t.z;
			element(r,3) = t.w;
		}
		
		void SetColumn(int c, const Vec4<T> & t) {
			element(0,c) = t.x;
			element(1,c) = t.y;
			element(2,c) = t.z;
			element(3,c) = t.w;
		}
		
		
		void GetRow(int r, Vec4<T> & t) const {
			t.x = element(r,0);
			t.y = element(r,1);
			t.z = element(r,2);
			t.w = element(r,3);
		}
		
		Vec4<T> GetRow(int r) const {
			Vec4<T> v; GetRow(r, v);
			return v;
		}
		
		void GetColumn(int c, Vec4<T> & t) const {
			t.x = element(0,c);
			t.y = element(1,c);
			t.z = element(2,c);
			t.w = element(3,c);
		}
		
		Vec4<T> GetColumn(int c) const {
			Vec4<T> v; GetColumn(c, v);
			return v;
		}
		
		Matrix4 Inverse() const	{
			Matrix4 minv;
			
			T r1[8], r2[8], r3[8], r4[8];
			T *s[4], *tmprow;
			
			s[0] = &r1[0];
			s[1] = &r2[0];
			s[2] = &r3[0];
			s[3] = &r4[0];
			
			register int i,j,p,jj;
			for(i=0;i<4;i++)
			{
				for(j=0;j<4;j++)
				{
					s[i][j] = element(i,j);
					if(i==j) s[i][j+4] = 1.0;
					else     s[i][j+4] = 0.0;
				}
			}
			T scp[4];
			for(i=0;i<4;i++)
			{
				scp[i] = T(fabs(s[i][0]));
				for(j=1;j<4;j++)
					if(T(fabs(s[i][j])) > scp[i]) scp[i] = T(fabs(s[i][j]));
				if(scp[i] == 0.0) return minv; // singular matrix!
			}
			
			int pivot_to;
			T scp_max;
			for(i=0;i<4;i++)
			{
				// select pivot row
				pivot_to = i;
				scp_max = T(fabs(s[i][i]/scp[i]));
				// find out which row should be on top
				for(p=i+1;p<4;p++)
					if(T(fabs(s[p][i]/scp[p])) > scp_max)
					{ scp_max = T(fabs(s[p][i]/scp[p])); pivot_to = p; }
				// Pivot if necessary
				if(pivot_to != i)
				{
					tmprow = s[i];
					s[i] = s[pivot_to];
					s[pivot_to] = tmprow;
					T tmpscp;
					tmpscp = scp[i];
					scp[i] = scp[pivot_to];
					scp[pivot_to] = tmpscp;
				}
				
				T mji;
				// perform gaussian elimination
				for(j=i+1;j<4;j++)
				{
					mji = s[j][i]/s[i][i];
					s[j][i] = 0.0;
					for(jj=i+1;jj<8;jj++)
						s[j][jj] -= mji*s[i][jj];
				}
			}
			if(s[3][3] == 0.0) return minv; // singular matrix!
			
			//
			// Now we have an upper triangular matrix.
			//
			//  x x x x | y y y y
			//  0 x x x | y y y y 
			//  0 0 x x | y y y y
			//  0 0 0 x | y y y y
			//
			//  we'll back substitute to get the inverse
			//
			//  1 0 0 0 | z z z z
			//  0 1 0 0 | z z z z
			//  0 0 1 0 | z z z z
			//  0 0 0 1 | z z z z 
			//
			
			T mij;
			for(i=3;i>0;i--)
			{
				for(j=i-1;j > -1; j--)
				{
					mij = s[j][i]/s[i][i];
					for(jj=j+1;jj<8;jj++)
						s[j][jj] -= mij*s[i][jj];
				}
			}
			
			for(i=0;i<4;i++)
				for(j=0;j<4;j++)
					minv(i,j) = s[i][j+4] / s[i][i];
			
			return minv;
		}
		
		
		Matrix4 Transpose() const {
			Matrix4 mtrans;
			
			for(int i=0;i<4;i++)
				for(int j=0;j<4;j++)
					mtrans(i,j) = element(j,i);		
			return mtrans;
		}
		
		Matrix4 & MultRight( const Matrix4 & b ) {
			Matrix4 mt(*this);
			SetValue(T(0));
			
			for(int i=0; i < 4; i++)
				for(int j=0; j < 4; j++)
					for(int c=0; c < 4; c++)
						element(i,j) += mt(i,c) * b(c,j);
			return *this;
		}    
		
		Matrix4 & MultLeft( const Matrix4 & b ) {
			Matrix4 mt(*this);
			SetValue(T(0));
			
			for(int i=0; i < 4; i++)
				for(int j=0; j < 4; j++)
					for(int c=0; c < 4; c++)
						element(i,j) += b(i,c) * mt(c,j);
			return *this;
		}
		
		
		// dst = M * src
		void MultMatrixVec( const Vec3<T> &src, Vec3<T> &dst ) const {
			T w = (
					  src.x * element(3,0) +
					  src.y * element(3,1) + 
					  src.z * element(3,2) +
					  element(3,3)          );
			
			assert(w != R3_ZERO);
			
			dst.x  = (
						 src.x * element(0,0) +
						 src.y * element(0,1) +
						 src.z * element(0,2) +
						 element(0,3)          ) / w;
			dst.y  = (
						 src.x * element(1,0) +
						 src.y * element(1,1) +
						 src.z * element(1,2) +
						 element(1,3)          ) / w;
			dst.z  = (
						 src.x * element(2,0) +
						 src.y * element(2,1) + 
						 src.z * element(2,2) +
						 element(2,3)          ) / w;
		}
		
		void MultMatrixVec( Vec3<T> & src_and_dst) const {
			MultMatrixVec(Vec3<T>(src_and_dst), src_and_dst); 
		}
		
		
		// dst = src * M
		void MultVecMatrix( const Vec3<T> &src, Vec3<T> &dst ) const {
			T w = (
					  src.x * element(0,3) +
					  src.y * element(1,3) +
					  src.z * element(2,3) +
					  element(3,3)          );
			
			assert(w != R3_ZERO);
			
			dst.x  = (
						 src.x * element(0,0) +
						 src.y * element(1,0) + 
						 src.z * element(2,0) + 
						 element(3,0)          ) / w;
			dst.y  = (
						 src.x * element(0,1) +
						 src.y * element(1,1) +
						 src.z * element(2,1) +
						 element(3,1)          ) / w;
			dst.z  = (
						 src.x * element(0,2) +
						 src.y * element(1,2) +
						 src.z * element(2,2) +
						 element(3,2)          ) / w;
		}
        
		
		void MultVecMatrix( Vec3<T> & src_and_dst) const {
			MultVecMatrix(Vec3<T>(src_and_dst), src_and_dst); 
		}
		
		// dst = M * src
		void MultMatrixVec( const Vec4<T> &src, Vec4<T> &dst ) const {
			dst.x  = (
						 src.x * element(0,0) +
						 src.y * element(0,1) +
						 src.z * element(0,2) +
						 src.w * element(0,3));
			dst.y  = (
						 src.x * element(1,0) +
						 src.y * element(1,1) +
						 src.z * element(1,2) +
						 src.w * element(1,3));
			dst.z  = (
						 src.x * element(2,0) +
						 src.y * element(2,1) + 
						 src.z * element(2,2) +
						 src.w * element(2,3));
			dst.w = (
						src.x * element(3,0) +
						src.y * element(3,1) + 
						src.z * element(3,2) +
						src.w * element(3,3));
		}
		
		void MultMatrixVec( Vec4<T> & src_and_dst) const {
			MultMatrixVec(Vec4<T>(src_and_dst), src_and_dst); 
		}
		
		
		// dst = src * M
		void MultVecMatrix( const Vec4<T> &src, Vec4<T> &dst ) const {
			dst.x  = (
						 src.x * element(0,0) +
						 src.y * element(1,0) + 
						 src.z * element(2,0) + 
						 src.w * element(3,0));
			dst.y  = (
						 src.x * element(0,1) +
						 src.y * element(1,1) +
						 src.z * element(2,1) +
						 src.w * element(3,1));
			dst.z  = (
						 src.x * element(0,2) +
						 src.y * element(1,2) +
						 src.z * element(2,2) +
						 src.w * element(3,2));
			dst.w = (
						src.x * element(0,3) +
						src.y * element(1,3) +
						src.z * element(2,3) +
						src.w * element(3,3));
		}
        
		
		void MultVecMatrix( Vec4<T> & src_and_dst) const {
			MultVecMatrix(Vec4<T>(src_and_dst), src_and_dst); 
		}
		
		
		// dst = M * src
		void MultMatrixDir( const Vec3<T> &src, Vec3<T> &dst ) const {
			dst.x  = (
						 src.x * element(0,0) +
						 src.y * element(0,1) +
						 src.z * element(0,2) ) ;
			dst.y  = ( 
						 src.x * element(1,0) +
						 src.y * element(1,1) +
						 src.z * element(1,2) ) ;
			dst.z  = ( 
						 src.x * element(2,0) +
						 src.y * element(2,1) + 
						 src.z * element(2,2) ) ;
		}
        
		
		void MultMatrixDir( Vec3<T> & src_and_dst) const {
			MultMatrixDir(Vec3<T>(src_and_dst), src_and_dst); 
		}
		
		
		// dst = src * M
		void MultDirMatrix( const Vec3<T> &src, Vec3<T> &dst ) const {
			dst.x  = ( 
						 src.x * element(0,0) +
						 src.y * element(1,0) +
						 src.z * element(2,0) ) ;
			dst.y  = ( 
						 src.x * element(0,1) +
						 src.y * element(1,1) +
						 src.z * element(2,1) ) ;
			dst.z  = (
						 src.x * element(0,2) +
						 src.y * element(1,2) + 
						 src.z * element(2,2) ) ;
		}
		
		
		void MultDirMatrix( Vec3<T> & src_and_dst) const {
			MultDirMatrix(Vec3<T>(src_and_dst), src_and_dst);
		}
		
		T * Ptr() {
			return m;
		}
		
		const T * Ptr() const {
			return m;
		}
		
		T & operator () (int row, int col) {
			return element(row,col); 
		}
		
		const T & operator () (int row, int col) const {
			return element(row,col); 
		}
		
		T & element (int row, int col) {
			return m[row | (col<<2)]; 
		}
		
		const T & element (int row, int col) const {
			return m[row | (col<<2)]; 
		}
		
		Matrix4 & operator *= ( const Matrix4 & mat ) {
			MultRight( mat );
			return *this;
		}
		
		Matrix4 & operator *= ( const T & r ) {
			for (int i = 0; i < 4; ++i)
			{
				element(0,i) *= r;
				element(1,i) *= r;
				element(2,i) *= r;
				element(3,i) *= r;
			}
			return *this;
		}
		
		Matrix4 & operator += ( const Matrix4 & mat ) {
			for (int i = 0; i < 4; ++i)
			{
				element(0,i) += mat.element(0,i);
				element(1,i) += mat.element(1,i);
				element(2,i) += mat.element(2,i);
				element(3,i) += mat.element(3,i);
			}
			return *this;
		}
		
		T m[16];
	};
	
	template <typename T> inline  
	Matrix4<T> operator * ( const Matrix4<T> & m1, const Matrix4<T> & m2 ) {
		static Matrix4<T> product;
		
		product = m1;
		product.MultRight(m2);
		
		return product;
	}
	
	template <typename T> inline
	bool operator ==( const Matrix4<T> &m1, const Matrix4<T> &m2 ) {
		return ( 
				m1(0,0) == m2(0,0) &&
				m1(0,1) == m2(0,1) &&
				m1(0,2) == m2(0,2) &&
				m1(0,3) == m2(0,3) &&
				m1(1,0) == m2(1,0) &&
				m1(1,1) == m2(1,1) &&
				m1(1,2) == m2(1,2) &&
				m1(1,3) == m2(1,3) &&
				m1(2,0) == m2(2,0) &&
				m1(2,1) == m2(2,1) &&
				m1(2,2) == m2(2,2) &&
				m1(2,3) == m2(2,3) &&
				m1(3,0) == m2(3,0) &&
				m1(3,1) == m2(3,1) &&
				m1(3,2) == m2(3,2) &&
				m1(3,3) == m2(3,3) );
	}
	
	template <typename T> inline
	bool operator != ( const Matrix4<T> & m1, const Matrix4<T> & m2 ) {
		return !( m1 == m2 ); 
	}  
	
	
	template <typename T> inline
	Vec3<T> operator * ( const Matrix4<T> & m, const Vec3<T> & v ) {
		Vec3<T> r;
		m.MultMatrixVec(v, r);
		return r;
	}
	
	template <typename T> inline
	Vec4<T> operator * ( const Matrix4<T> & m, const Vec4<T> & v ) {
		Vec4<T> r;
		m.MultMatrixVec(v, r);
		return r;
	}
	
	
	template <typename T> inline
	Vec3<T> operator * ( const Vec3<T> & v, const Matrix4<T> & m ) {
		Vec3<T> r;
		m.MultVecMatrix(v, r);
		return r;
	}
	
	template <typename T> inline
	Vec4<T> operator * ( const Vec4<T> & v, const Matrix4<T> & m ) {
		Vec4<T> r;
		m.MultVecMatrix(v, r);
		return r;
	}
	
	
	template <typename T> inline
	Matrix4<T> ToMatrix4( const Matrix3<T> & m3 ) {
		Matrix4<T> m4;
		m4( 0, 0 ) = m3( 0, 0 );   m4( 0, 1 ) = m3( 0, 1 );   m4( 0, 2 ) = m3( 0, 2 );   m4( 0, 3 ) = 0.0f;
		m4( 1, 0 ) = m3( 1, 0 );   m4( 1, 1 ) = m3( 1, 1 );   m4( 1, 2 ) = m3( 1, 2 );   m4( 1, 3 ) = 0.0f;
		m4( 2, 0 ) = m3( 2, 0 );   m4( 2, 1 ) = m3( 2, 1 );   m4( 2, 2 ) = m3( 2, 2 );   m4( 2, 3 ) = 0.0f;
		m4( 3, 0 ) = 0.0f;         m4( 3, 1 ) = 0.0f;         m4( 3, 2 ) = 0.0f;         m4( 3, 3 ) = 1.0f;
		return m4;
	}
	
	
	template <typename T> inline
	Matrix3<T> ToMatrix3( const Matrix4<T> & m4 ) {
		Matrix3<T> m3;
		m3( 0, 0 ) = m4( 0, 0 );   m3( 0, 1 ) = m4( 0, 1 );   m3( 0, 2 ) = m4( 0, 2 );
		m3( 1, 0 ) = m4( 1, 0 );   m3( 1, 1 ) = m4( 1, 1 );   m3( 1, 2 ) = m4( 1, 2 );
		m3( 2, 0 ) = m4( 2, 0 );   m3( 2, 1 ) = m4( 2, 1 );   m3( 2, 2 ) = m4( 2, 2 );
		return m3;
	}
	
	
	
	
	
	
	
	template < typename T>
	class Quaternion {
	public:
		
		Quaternion() {
			q[0]=q[1]=q[2]=0.0;
			q[3]=R3_ONE;
		}
		
		Quaternion( const T v[4] ) {
			SetValue( v );
		}
		
		
		Quaternion( T q0, T q1, T q2, T q3 ) {
			SetValue( q0, q1, q2, q3 );
		}
		
		
		Quaternion( const Matrix4<T> & m )	{
			SetValue( m );
		}
		
		
		Quaternion( const Vec3<T> &axis, T radians ) {
			SetValue( axis, radians );
		}
		
		
		Quaternion( const Vec3<T> &rotateFrom, const Vec3<T> &rotateTo ) {
			SetValue( rotateFrom, rotateTo );
		}
		
		Quaternion( const Vec3<T> & fromLook, const Vec3<T> & fromUp,
				   const Vec3<T> & toLook, const Vec3<T>& toUp) {
			SetValue(fromLook, fromUp, toLook, toUp);
		}
		
		const T * GetValue() const {
			return  &q[0];
		}
		
		void GetValue( T &q0, T &q1, T &q2, T &q3 ) const {
			q0 = q[0];
			q1 = q[1];
			q2 = q[2];
			q3 = q[3];
		}
		
		Quaternion & SetValue( T q0, T q1, T q2, T q3 ) {
			q[0] = q0;
			q[1] = q1;
			q[2] = q2;
			q[3] = q3;
			counter = 0;
			return *this;
		}
		
		void GetValue( Vec3<T> &axis, T &radians ) const {
			radians = T(acos( q[3] ) * R3_TWO);
			if ( radians == R3_ZERO )
				axis = Vec3<T>( 0.0, 0.0, 1.0 );
			else
			{
				axis.x = q[0];
				axis.y = q[1];
				axis.z = q[2];
				axis.Normalize();
			}
		}
		
		void GetValue( Matrix3<T> & m ) const {
			T s, xs, ys, zs, wx, wy, wz, xx, xy, xz, yy, yz, zz;
			
			T norm = q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3];
			
			s = (r3_equivalent(norm,R3_ZERO)) ? R3_ZERO : ( R3_TWO / norm );
			
			xs = q[0] * s;
			ys = q[1] * s;
			zs = q[2] * s;
			
			wx = q[3] * xs;
			wy = q[3] * ys;
			wz = q[3] * zs;
			
			xx = q[0] * xs;
			xy = q[0] * ys;
			xz = q[0] * zs;
			
			yy = q[1] * ys;
			yz = q[1] * zs;
			zz = q[2] * zs;
			
			m(0,0) = T( R3_ONE - ( yy + zz ));
			m(1,0) = T ( xy + wz );
			m(2,0) = T ( xz - wy );
			
			m(0,1) = T ( xy - wz );
			m(1,1) = T ( R3_ONE - ( xx + zz ));
			m(2,1) = T ( yz + wx );
			
			m(0,2) = T ( xz + wy );
			m(1,2) = T ( yz - wx );
			m(2,2) = T ( R3_ONE - ( xx + yy ));
		}
		
		void GetValue( Matrix4<T> & m ) const {
			Matrix3<T> m3;
			GetValue( m3 );
			m = ToMatrix4( m3 );
		}
		
		Matrix3<T> GetMatrix3() const {
			Matrix3<T> m;
			GetValue( m );
			return m;
		}

		Matrix4<T> GetMatrix4() const {
			Matrix4<T> m;
			GetValue( m );
			return m;
		}
		
		Quaternion & SetValue( const T * qp ) {
			memcpy(q,qp,sizeof(T) * 4);
			
			counter = 0;
			return *this;
		}
		
		Quaternion & SetValue( const Matrix4<T> & m ) {
			T tr, s;
			int i, j, k;
			const int nxt[3] = { 1, 2, 0 };
			
			tr = m(0,0) + m(1,1) + m(2,2);
			
			if ( tr > R3_ZERO )
			{
				s = T(sqrt( tr + m(3,3) ));
				q[3] = T ( s * 0.5 );
				s = T(0.5) / s;
				
				q[0] = T ( ( m(1,2) - m(2,1) ) * s );
				q[1] = T ( ( m(2,0) - m(0,2) ) * s );
				q[2] = T ( ( m(0,1) - m(1,0) ) * s );
			}
			else
			{
				i = 0;
				if ( m(1,1) > m(0,0) )
					i = 1;
				
				if ( m(2,2) > m(i,i) )
					i = 2;
				
				j = nxt[i];
				k = nxt[j];
				
				s = T(sqrt( ( m(i,j) - ( m(j,j) + m(k,k) )) + R3_ONE ));
				
				q[i] = T ( s * 0.5 );
				s = T(0.5 / s);
				
				q[3] = T ( ( m(j,k) - m(k,j) ) * s );
				q[j] = T ( ( m(i,j) + m(j,i) ) * s );
				q[k] = T ( ( m(i,k) + m(k,i) ) * s );
			}
			
			counter = 0;
			return *this;
		}
		
		Quaternion & SetValue( const Vec3<T> &axis, T theta ) {
			T sqnorm = axis.SquareNorm();
			
			if (sqnorm <= R3_EPSILON)
			{
				// axis too small.
				x = y = z = 0.0;
				w = R3_ONE;
			} 
			else 
			{
				theta *= T(0.5);
				T sin_theta = T(sin(theta));
				
				if (!r3_equivalent(sqnorm,R3_ONE)) 
					sin_theta /= T(sqrt(sqnorm));
				x = sin_theta * axis.x;
				y = sin_theta * axis.y;
				z = sin_theta * axis.z;
				w = T(cos(theta));
			}
			return *this;
		}
		
		Quaternion & SetValue( const Vec3<T> & rotateFrom, const Vec3<T> & rotateTo ) {
			Vec3<T> p1, p2;
			T alpha;
			
			p1 = rotateFrom; 
			p1.Normalize();
			p2 = rotateTo;  
			p2.Normalize();
			
			alpha = p1.Dot(p2);
			
			if(r3_greater_than(alpha,R3_ONE))
			{ 
				*this = Identity(); 
				return *this; 
			}
			
			// ensures that the anti-parallel case leads to a positive dot
			if(r3_less_than(alpha,-R3_ONE))
			{
				Vec3<T> v ;
				
				if( p1.x != p1.y )
					v = Vec3<T>(p1.y, p1.x, p1.z);
				else
					v = Vec3<T>(p1.z, p1.y, p1.x);
				
				v -= p1 * p1.Dot(v);
				v.Normalize();
				
				SetValue(v, R3_PI);
				return *this;
			}
			
			p1 = p1.Cross(p2);  
			p1.Normalize();
			SetValue(p1,T(acos(alpha)));
			
			counter = 0;
			return *this;
		}
		
		Quaternion & SetValue( const Vec3<T> & fromLook, const Vec3<T> & fromUp,
							   const Vec3<T> & toLook, const Vec3<T> & toUp) {
			Quaternion rLook = Quaternion(fromLook, toLook);
			
			Vec3<T> rotatedFromUp(fromUp);
			rLook.MultVec(rotatedFromUp);
			
			Quaternion rTwist = Quaternion(rotatedFromUp, toUp);
			
			*this = rTwist;
			*this *= rLook;
			return *this;
		}
		
		Quaternion & operator *= ( const Quaternion & qr ) {
			Quaternion ql(*this);
			
			w = ql.w * qr.w - ql.x * qr.x - ql.y * qr.y - ql.z * qr.z;
			x = ql.w * qr.x + ql.x * qr.w + ql.y * qr.z - ql.z * qr.y;
			y = ql.w * qr.y + ql.y * qr.w + ql.z * qr.x - ql.x * qr.z;
			z = ql.w * qr.z + ql.z * qr.w + ql.x * qr.y - ql.y * qr.x;
			
			counter += qr.counter;
			counter++;
			CounterNormalize();
			return *this;
		}
		
		void Normalize() {
			T rnorm = R3_ONE / T(sqrt(w * w + x * x + y * y + z * z));
			if (r3_equivalent(rnorm, R3_ZERO))
				return;
			x *= rnorm;
			y *= rnorm;
			z *= rnorm;
			w *= rnorm;
			counter = 0;
		}
		
	
		bool Equals( const Quaternion & r, T tolerance ) const {
			T t;
			
			t = (
				 (q[0]-r.q[0])*(q[0]-r.q[0]) +
				 (q[1]-r.q[1])*(q[1]-r.q[1]) +
				 (q[2]-r.q[2])*(q[2]-r.q[2]) +
				 (q[3]-r.q[3])*(q[3]-r.q[3]) );
			if(t > R3_EPSILON) 
				return false;
			return 1;
		}
		
		Quaternion & Conjugate() {
			q[0] *= -R3_ONE;
			q[1] *= -R3_ONE;
			q[2] *= -R3_ONE;
			return *this;
		}
		
		Quaternion & Invert() {
			return Conjugate();
		}
		
		Quaternion Inverse() const {
			Quaternion r = *this;
			return r.Invert();
		}
		
		//
		// Quaternion multiplication with cartesian vector
		// v' = q*v*q(star)
		//
		void MultVec( const Vec3<T> &src, Vec3<T> &dst ) const {
			T v_coef = w * w - x * x - y * y - z * z;                     
			T u_coef = R3_TWO * (src.x * x + src.y * y + src.z * z);  
			T c_coef = R3_TWO * w;                                       
			
			dst.x = v_coef * src.x + u_coef * x + c_coef * (y * src.z - z * src.y);
			dst.y = v_coef * src.y + u_coef * y + c_coef * (z * src.x - x * src.z);
			dst.z = v_coef * src.z + u_coef * z + c_coef * (x * src.y - y * src.x);
		}
		
		void MultVec( Vec3<T> & src_and_dst) const {
			MultVec(Vec3<T>(src_and_dst), src_and_dst);
		}
		
		void ScaleAngle( T scaleFactor ) {
			Vec3<T> axis;
			T radians;
			
			GetValue(axis, radians);
			radians *= scaleFactor;
			SetValue(axis, radians);
		}
		
		static Quaternion Slerp( const Quaternion & p, const Quaternion & q, T alpha ) {
			Quaternion r;
			
			T cos_omega = p.x * q.x + p.y * q.y + p.z * q.z + p.w * q.w;
			// if B is on opposite hemisphere from A, use -B instead
			
			int bflip;
			if ( ( bflip = (cos_omega < R3_ZERO)) )
				cos_omega = -cos_omega;
			
			// complementary interpolation parameter
			T beta = R3_ONE - alpha;     
			
			if(cos_omega <= R3_ONE - R3_EPSILON)
				return p;
			
			T omega = T(acos(cos_omega));
			T one_over_sin_omega = R3_ONE / T(sin(omega));
			
			beta    = T(sin(omega*beta)  * one_over_sin_omega);
			alpha   = T(sin(omega*alpha) * one_over_sin_omega);
			
			if (bflip)
				alpha = -alpha;
			
			r.x = beta * p.q[0]+ alpha * q.q[0];
			r.y = beta * p.q[1]+ alpha * q.q[1];
			r.z = beta * p.q[2]+ alpha * q.q[2];
			r.w = beta * p.q[3]+ alpha * q.q[3];
			return r;
		}
		
		static Quaternion Identity() {
			static Quaternion ident( Vec3<T>( 0.0, 0.0, 0.0 ), R3_ONE );
			return ident;
		}
		
		T & operator []( int i ) {
			assert(i < 4);
			return q[i];
		}
		
		const T & operator []( int i ) const {
			assert(i < 4);
			return q[i];
		}
		
		union {
			struct { T q[4]; };
			struct { T x, y, z, w; };
		};
		
	protected:
		
		void CounterNormalize() {
			if (counter > R3_QUATERNION_NORMALIZATION_THRESHOLD)
				Normalize();
		}
		
		// renormalization counter
		unsigned char counter;
	};
	
	template <typename T> inline
	bool operator == ( const Quaternion<T> & q1, const Quaternion<T> & q2 ) {
		return (r3_equivalent(q1.x, q2.x) &&
				r3_equivalent(q1.y, q2.y) &&
				r3_equivalent(q1.z, q2.z) &&
				r3_equivalent(q1.w, q2.w) );
	}
	
	template <typename T> inline
	bool operator != ( const Quaternion<T> & q1, const Quaternion<T> & q2 ) { 
		return ! ( q1 == q2 ); 
	}
	
	template <typename T> inline
	Quaternion<T> operator * ( const Quaternion<T> & q1, const Quaternion<T> & q2 ) {	
		Quaternion<T> r(q1); 
		r *= q2; 
		return r; 
	}
	
	
	
	
	
	
	
	
	
	template < typename T>
	class Plane {
	public:
		
		Plane()	{
			planedistance = 0.0;
			planenormal.SetValue( 0.0, 0.0, 1.0 );
		}
		
		
		Plane( const Vec3<T> &p0, const Vec3<T> &p1, const Vec3<T> &p2 ) {
			Vec3<T> v0 = p1 - p0;
			Vec3<T> v1 = p2 - p0;
			planenormal = v0.Cross(v1);  
			planenormal.Normalize();
			planedistance = p0.Dot(planenormal);
		}
		
		Plane( const Vec3<T> &normal, T distance ) {
			planedistance = distance;
			planenormal = normal;
			planenormal.Normalize();
		}
		
		Plane( const Vec3<T> &normal, const Vec3<T> &point ) {
			planenormal = normal;
			planenormal.Normalize();
			planedistance = point.Dot(planenormal);
		}
		
		void Offset( T d ) {
			planedistance += d;
		}
		
		bool Intersect( const Line<T> &l, Vec3<T> &intersection ) const {
			Vec3<T> pos, dir;
			Vec3<T> pn = planenormal;
			T pd = planedistance;
			
			pos = l.GetPosition();
			dir = l.GetDirection();
			
			if(dir.Dot(pn) == 0.0) return 0;
			pos -= pn*pd;
			// now we're talking about a Plane passing through the origin
			if(pos.Dot(pn) < 0.0) pn.Negate();
			if(dir.Dot(pn) > 0.0) dir.Negate();
			Vec3<T> ppos = pn * pos.Dot(pn);
			pos = (ppos.Length()/dir.Dot(-pn))*dir;
			intersection = l.GetPosition();
			intersection += pos;
			return 1;
		}
		
		void Transform( const Matrix4<T> &matrix ) {
			Matrix4<T> invtr = matrix.Inverse();
			invtr = invtr.Transpose();
			
			Vec3<T> pntOnPlane = planenormal * planedistance;
			Vec3<T> newPntOnPlane;
			Vec3<T> newnormal;
			
			invtr.MultDirMatrix(planenormal, newnormal);
			matrix.MultVecMatrix(pntOnPlane, newPntOnPlane);
			
			newnormal.Normalize();
			planenormal = newnormal;
			planedistance = newPntOnPlane.Dot(planenormal);
		}
		
		bool IsInHalfSpace( const Vec3<T> &point ) const {
			
			if(( point.Dot(planenormal) - planedistance) < 0.0)
				return 0;
			return 1;
		}
		
		
		T Distance( const Vec3<T> & point ) const {
			return planenormal.Dot(point - planenormal*planedistance);
		}
		
		const Vec3<T> &GetNormal() const {
			return planenormal;
		}
		
		
		T GetDistanceFromOrigin() const {
			return planedistance;
		}
		
		
		
		//protected:
		Vec3<T> planenormal;
		T planedistance;
	};
	
	template <typename T> inline
	bool operator == (const Plane<T> & p1, const Plane<T> & p2 ) {
		return (  p1.planedistance == p2.planedistance && p1.planenormal == p2.planenormal);
	}
	
	template <typename T> inline
	bool operator != ( const Plane<T> & p1, const Plane<T> & p2 ) {
		return  ! (p1 == p2);
	}
	
	
	
	
	
	
	// some useful constructors / functions
	
    // inverse of camera_lookat
	template <typename T> inline
    Matrix4<T> ObjectLookAt(const Vec3<T> & from, const Vec3<T> & to, const Vec3<T> & Up)
    {
        Vec3<T> look = to - from;
        look.Normalize();
        Vec3<T> up(Up);
        up -= look * look.Dot(up);
        up.Normalize();
		
        Quaternion<T> r( Vec3<T>(0,0,-1), Vec3<T>(0,1,0), look, up);
        Matrix4<T> m;
        r.GetValue(m);
        m.SetTranslate(from);
        return m;
    }
	
    // inverse of object_lookat
	template <typename T> inline
    Matrix4<T> CameraLookAt(const Vec3<T> & eye, const Vec3<T> & lookpoint, const Vec3<T> & Up)
    {
        Vec3<T> look = lookpoint - eye;
        look.Normalize();
        Vec3<T> up(Up);
        up -= look * look.Dot(up);
        up.Normalize();
		
        Matrix4<T> t;
        t.SetTranslate(-eye);
		
        Quaternion<T> r( Vec3<T>(0,0,-1), Vec3<T>(0,1,0), look, up);
        r.Invert();
        Matrix4<T> rm;
        r.GetValue(rm);
        return rm * t;
    }
	
	
	
	template <typename T> inline
	Matrix4<T> Frustum( T left, T right, T bottom, T top, T zNear, T zFar )
    {
        Matrix4<T> m;
        m.MakeIdentity();
		
        m(0,0) = (2*zNear) / (right - left);
        m(0,2) = (right + left) / (right - left);
		
        m(1,1) = (2*zNear) / (top - bottom);
        m(1,2) = (top + bottom) / (top - bottom);
		
        m(2,2) = -(zFar + zNear) / (zFar - zNear);
        m(2,3) = -2*zFar*zNear / (zFar - zNear);
		
        m(3,2) = -1;
        m(3,3) = 0;
		
        return m;
    }
	
	template <typename T> inline
	Matrix4<T> FrustumInverse( T left, T right,	T bottom, T top, T zNear, T zFar)
    {
        Matrix4<T> m;
        m.MakeIdentity();
		
        m(0,0) = (right - left) / (2 * zNear);
        m(0,3) = (right + left) / (2 * zNear);
		
        m(1,1) = (top - bottom) / (2 * zNear);
        m(1,3) = (top + bottom) / (2 * zNear);
		
        m(2,2) = 0;
        m(2,3) = -1;
		
        m(3,2) = -(zFar - zNear) / (2 * zFar * zNear);
        m(3,3) = (zFar + zNear) / (2 * zFar * zNear);
		
        return m;
    }
	
	template <typename T> inline
    Matrix4<T> Perspective( T fovy, T aspect, T zNear, T zFar)
    {
        T tangent = (T)tan(ToRadians(fovy/T(2.0)));
        T y = tangent * zNear;
        T x = aspect * y;
        return Frustum(-x, x, -y, y, zNear, zFar);
    }
	
	template <typename T> inline
    Matrix4<T> PerspectiveInverse( T fovy, T aspect, T zNear, T zFar)
    {
        T tangent = (T)tan(ToRadians(fovy/T(2.0)));
        T y = tangent * zNear;
        T x = aspect * y;
        return FrustumInverse(-x, x, -y, y, zNear, zFar);
    }
	
	template <typename T> inline 
	Matrix4<T> Ortho( T left, T right, T bottom, T top, T zNear, T zFar ) {
		Matrix4<T> m;
		Vec3<T> s( 1 / ( right - left ), 1 / ( top - bottom ), 1 / ( zFar - zNear ) );
		m.SetScale( s * T(2) );
		m.SetTranslate( s * Vec3<T>( -( right + left ), -( top + bottom ), zFar + zNear ) );
		return m;
	}
		
	template <typename T> inline 
	Matrix4<T> OrthoInverse( T left, T right, T bottom, T top, T zNear, T zFar ) {
		// could be done with a formula, but I'm being lazy
		Matrix4<T> m = Ortho( left, right, bottom, top, zNear, zFar );
		return m.Inverse();
	}
	
	
	
	// make common typedefs...
	typedef Vec2<int> Vec2i;
	typedef Vec2<float> Vec2f;
	typedef Vec2<double> Vec2d;
	typedef Vec3<int> Vec3i;
	typedef Vec3<float> Vec3f;
	typedef Vec3<double> Vec3d;
	typedef Vec4<float>  Vec4f;
	typedef Vec4<double> Vec4d;
	typedef Vec4<int> Vec4i;
	typedef Vec4< unsigned char > Vec4ub;
	typedef Quaternion<float>  Quaternionf;
	typedef Quaternion<double>  Quaterniond;
	typedef Quaternion<float>  Rotationf;
	typedef Quaternion<double>  Rotationd;
	typedef Line<float>  Linef;
	typedef Line<double>  Lined;
	typedef LineSegment2<float>  LineSegment2f;
	typedef LineSegment2<double>  LineSegment2d;
	typedef Plane<float>  Planef;
	typedef Plane<double>  Planed;
	typedef Matrix3<float> Matrix3f;
	typedef Matrix3<double> Matrix3d;
	typedef Matrix3<int> Matrix3i;	
	typedef Matrix4<float>  Matrix4f;
	typedef Matrix4<double>  Matrix4d;
	
}  // namespace r3



#endif


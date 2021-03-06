/*
Pup routines for STL classes.

After including this header, you can parameter-marshall
an variable consisting of STL vectors, lists, maps,
strings, or pairs.

This includes variables of type "std::list<int>", or even
"std::map<double, std::vector<std::string> >".

NOT included are the rarer types like valarray or slice, 
vector<bool>, set or multiset, or deque.

Orion Sky Lawlor, uiuc.edu@acm.org, 7/22/2002
*/
#ifndef _UIUC_CHARM_PUP_STL_H
#define _UIUC_CHARM_PUP_STL_H

/*It's kind of annoying that we have to drag all these headers in
  just so the std:: parameter declarations will compile.
 */
#include <vector>
#include <list>
#include <map>
#include <string>
#include <utility> /*for std::pair*/
#include "pup.h"

/*************** Simple classes ***************/

template <class A,class B> 
inline void PUPt(PUP::er &p,typename std::pair<const A,B> &v)
{
  p.syncComment(PUP::sync_index);
  PUPt(p,*(A *)&v.first); /* cast away constness on A */
  p.syncComment(PUP::sync_item);
  PUPt(p,v.second);
}
template <class charType> 
inline void PUPt(PUP::er &p,typename std::basic_string<charType> &v)
{
  int nChar=v.length();
  PUPt(p,nChar);
  if (p.isUnpacking()) { //Unpack to temporary buffer
    charType *buf=new charType[nChar];
    p(buf,nChar);
    v=std::basic_string<charType>(buf,nChar);
    delete[] buf;
  }
  else /*packing*/ { //Do packing in-place from data
    //Have to cast away constness here
    p((charType *)v.data(),nChar);
  }
}
inline void PUPt(PUP::er &p,std::string &v)
{
  p.syncComment(PUP::sync_begin_object,"std::string");
  int nChar=v.length();
  PUPt(p,nChar);
  if (p.isUnpacking()) { //Unpack to temporary buffer
    char *buf=new char[nChar];
    p(buf,nChar);
    v=std::basic_string<char>(buf,nChar);
    delete[] buf;
  }
  else /*packing*/ { //Do packing in-place from data
    //Have to cast away constness here
    p((char *)v.data(),nChar);
  }
  p.syncComment(PUP::sync_end_object);
}

/**************** Containers *****************/

//Impl. util: pup the length of a container
template <class container>
inline int PUP_stl_container_size(PUP::er &p,container &c) {
  int nElem=c.size();
  PUPt(p,nElem);
  return nElem; 
}

//Impl. util: pup each current item of a container (no allocation)
template <class container>
inline void PUP_stl_container_items(PUP::er &p,container &c) {
  for (typename container::iterator it=c.begin();
       it!=c.end();
       ++it) {
    p.syncComment(PUP::sync_item);
    PUPt(p,*it);  
  }
}

template <class container,class dtype>
inline void PUP_stl_container(PUP::er &p,container &c,dtype *data_type) {
  p.syncComment(PUP::sync_begin_array);
  int nElem=PUP_stl_container_size(p,c);
  if (p.isUnpacking()) 
  { //Unpacking: Extract each element and push_back:
    c.resize(0);
    for (int i=0;i<nElem;i++) {
      p.syncComment(PUP::sync_item);
      dtype n;
      PUPt(p,n);
      c.push_back(n);
    } 
  }
  else PUP_stl_container_items(p,c);
  p.syncComment(PUP::sync_end_array);
}

//Map objects don't have a "push_back", while vector and list
//  don't have an "insert", so PUP_stl_map isn't PUP_stl_container
template <class container,class dtype>
inline void PUP_stl_map(PUP::er &p,container &c,dtype *data_type) {
  p.syncComment(PUP::sync_begin_list);
  int nElem=PUP_stl_container_size(p,c);
  if (p.isUnpacking()) 
  { //Unpacking: Extract each element and insert:
    for (int i=0;i<nElem;i++) {
      dtype n;
      PUPt(p,n);
      c.insert(n);
    } 
  }
  else PUP_stl_container_items(p,c);
  p.syncComment(PUP::sync_end_list);
}

template <class T> 
inline void PUPt(PUP::er &p,typename std::vector<T> &v)
  { PUP_stl_container(p,v,(T *)0); }
template <class T> 
inline void PUPt(PUP::er &p,typename std::list<T> &v)
  { PUP_stl_container(p,v,(T *)0); }

template <class V,class T,class Cmp> 
inline void PUPt(PUP::er &p,typename std::map<V,T,Cmp> &m)
  { PUP_stl_map(p,m,(std::pair<const V,T> *)0); }
template <class V,class T,class Cmp> 
inline void PUPt(PUP::er &p,typename std::multimap<V,T,Cmp> &m)
  { PUP_stl_map(p,m,(std::pair<const V,T> *)0); }


#endif

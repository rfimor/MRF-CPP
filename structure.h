#ifndef STRUCTURE_H
#define STRUCTURE_H

#include <queue>
#include <stdio.h>

template <class T>
class Node
{

  public:
    T data;
    Node<T> *next;
    Node<T> *prev;

    //constructor
    Node(const T &item, Node<T> *ptrnext = NULL,
		 Node<T> *ptrprev = NULL)
		 : data(item), next(ptrnext) , prev(ptrprev) {}
    ~Node(void) {}

} ;

template <class T>
Node<T> * GetNode(const T &item, Node<T> *nextPtr = NULL,
                 Node<T> *prevPtr = NULL)
{
      Node<T> *newNode;
      newNode = new Node<T>(item,nextPtr,prevPtr);
      return newNode;
}
                 
template <class T>
class List
{
  public:
    Node<T> *head;
    Node<T> *tail;
  private:
    long int length;
  
  public:

    //constructor, construct an empty list
    List();

    virtual ~List(void);

    void InsertFront(const T &item);
    Node<T>* ListSearch(const T &item);
    void ListRemove(Node<T> * itemPtr);
    void ListDelete(Node<T> * itemPtr);
    void InsertRear(const T &item);
    void InsertBefore(const T &item, Node<T> *itemPtr);
    void InsertAfter(const T &item, Node<T> *itemPtr);
    long int ListLength();
    
    T* toArray(void);
} ;

template <class T>
class ListPtr
{
  public:
    Node<T*> *head;
    Node<T*> *tail;
  private:
    long int length;

  public:
    //constructor, construct an empty list
    ListPtr();

    virtual ~ListPtr(void);

    void InsertFront(T * &item);
    Node<T*>* ListSearch(T *item);
    void ListDelete(Node<T*> *itemPtr);
    void InsertRear(T * &item);
    void InsertBefore(T* &item, Node<T*> *itemPtr);
    void InsertAfter(T* &item, Node<T*> *itemPtr);
    void ListEmpty(void);
    void DeleteData(void);
    long int ListLength();
    
    T** toArray(void);
} ;

template <class T>
class Queue
{
   private:
     List<T>* ListQ;

   public:

     Queue(void);
     ~Queue(void);
     void Enqueue(const T &item);
     bool Dequeue(T &rval);
     long int QueueLength(void);
} ;


template <class T>
class QueuePtr
{
   private:
     ListPtr<T>* ListQ;

   public:

     QueuePtr(void);
     ~QueuePtr(void);
     void Enqueue(T * &item);
     bool Dequeue(T * &rval);
     long int QueueLength(void);
     void DeleteData(void);
} ;

class Point3D
{
   public:

       double hfield; 
       double delta;
       double energy;

       Point3D(void);
       Point3D(double h, double d, double e);
       bool operator == (const Point3D & p);
       void operator = (const Point3D & p);
} ;

class Point2D
{
   public:

       double hfield;
       double delta;

       Point2D(void);
       Point2D(double h,double d);
       bool operator == (const Point2D & p);
       void operator = (const Point2D & p);
       void operator = (const Point3D & p);
} ;

class State
{
   public:
	
       double couple;
       double rndmag;
       double mag;

       State(void) {}
       State(double c,double r,double m);
       State(const State & p);

       bool operator == (const State & p);
       void operator = (const State & p);
       bool OnState(const Point3D *p,const double EPS);
       bool OnState(const Point3D *p, double &diff, const double EPS);
       double Energy(const Point2D* test);
} ;

template <class T> void List<T>::InsertBefore(const T &item, Node<T> *itemPtr) {
	if (itemPtr == head) InsertFront(item);
	else {
		itemPtr->prev = itemPtr->prev->next = GetNode(item, itemPtr, itemPtr->prev);
		length++;
	}
}

template <class T> void List<T>::InsertAfter(const T &item, Node<T> *itemPtr) {
	if (itemPtr == tail) InsertRear(item);
	else {
		itemPtr->next = itemPtr->next->prev = GetNode(item, itemPtr->next, itemPtr);
		length++;
	}
}

template <class T> void ListPtr<T>::InsertBefore(T* &item, Node<T*> *itemPtr) {
	if (itemPtr == head) InsertFront(item);
	else {
		itemPtr->prev = itemPtr->prev->next = GetNode(item, itemPtr, itemPtr->prev);
		length++;
	}
}

template <class T> void ListPtr<T>::InsertAfter(T* &item, Node<T*> *itemPtr) {
	if (itemPtr == tail) InsertRear(item);
	else {
		itemPtr->next = itemPtr->next->prev = GetNode(item, itemPtr->next, itemPtr);
		length++;
	}
}

template <class T> T* List<T>::toArray() {
	T * ret = new T[length];
	Node<T> *curr = head;
	for (int i=0; i<length; i++) {
		ret[i] = curr->data;
		curr = curr->next;
	}
	return ret;
}

template <class T> T** ListPtr<T>::toArray() {
	T ** ret = new T*[length];
	Node<T*> *curr = head;
	for (int i=0; i<length; i++) {
		ret[i] = curr->data;
		curr = curr->next;
	}
	return ret;
}

template <class T> List<T> :: List()
{
       head = NULL;
       tail = NULL;
       length = 0;
}

template <class T> long int List<T> :: ListLength()
{
	return length;
}

template <class T> long int ListPtr<T> :: ListLength()
{
	return length;
}

template <class T> ListPtr<T> :: ListPtr()
{
       head = NULL;
       tail = NULL;
       length = 0;
}

    //Insert item to the front of the list
template <class T> void List<T> :: InsertFront(const T &item)
{
    if (head == NULL) {
       tail = head = GetNode(item,head);
    }
    else {
       head = head->prev = GetNode(item,head);
    }
    length++;
}

    //search an item x in the list.
	//if found output is pointer of node that contains x
    // otherwise NULL
template <class T> Node<T>* List<T> :: ListSearch(const T &item)
{
     Node<T>* x;
     x = head;
     while (x != NULL && (x->data) != item) {
         x = x->next;
     }
     return x;
}

     //remove an node from the list
template <class T> void List<T> :: ListRemove(Node<T> * itemPtr)
{
     if (itemPtr->prev != NULL) {
        (itemPtr->prev)->next = itemPtr->next;
     }
     else {
        head = itemPtr->next;
     }

     if (itemPtr->next != NULL) {
        (itemPtr->next)->prev = itemPtr->prev;
     }
     else {
        tail = itemPtr->prev;
     }
     length--;
}

    //remove an node from the list, also delete the node
template <class T> void List<T> :: ListDelete(Node<T> * itemPtr)
{
     if (itemPtr->prev != NULL) {
        (itemPtr->prev)->next = itemPtr->next;
     }
     else {
        head = itemPtr->next;
     }

     if (itemPtr->next != NULL) {
        (itemPtr->next)->prev = itemPtr->prev;
     }
     else {
        tail = itemPtr->prev;
     }
     delete itemPtr;
     length--;
}

    //Insert item to the rear of the list
template <class T> void List<T> :: InsertRear(const T &item)
{
     Node<T> *tempNull;
     tempNull = NULL;

     if (head == NULL) {
        InsertFront(item);
     }
     else  {
        tail = tail->next = GetNode(item,tempNull,tail);
        length++;
     }
}

template <class T> List<T> :: ~List()
{
     while (head != NULL) {
           ListDelete(head);
     }
}
//------------------------------------------------------------------------

    //Insert item to the front of the list
template <class T> void ListPtr<T> :: InsertFront(T * &item)
{
    if (head == NULL) {
       tail = head = GetNode(item,head);
    }
    else {
       head = head->prev = GetNode(item,head);
    }
    length++;
}

    //search an item x in the list.
	//if found output is pointer of node that contains x
    // otherwise NULL
template <class T> Node<T*>* ListPtr<T> :: ListSearch(T *item)
{
     Node<T*>* x;
     x = head;
     while (x != NULL && (x->data) != item) {
         x = x->next;
     }
     return x;
}


    //remove an item from the list, NEVER delete the DATA in  nodes here!
template <class T> void ListPtr<T> :: ListDelete(Node<T*> *itemPtr)
{
     if (itemPtr->prev != NULL) {
        (itemPtr->prev)->next = itemPtr->next;
     }
     else {
        head = itemPtr->next;
     }

     if (itemPtr->next != NULL) {
        (itemPtr->next)->prev = itemPtr->prev;
     }
     else {
        tail = itemPtr->prev;
     }

     delete itemPtr;
     length--;
}

    //Insert item to the rear of the list
template <class T> void ListPtr<T> :: InsertRear(T * &item)
{
     if (head == NULL) {
        InsertFront(item);
        return;
     }
     else  {
        tail = tail->next = GetNode(item,(Node<T*> *)NULL,tail);
         length++;
     }
}

template <class T> ListPtr<T> :: ~ListPtr()
{
     while (head != NULL) {
           ListDelete(head);
     }
}

template <class T> void ListPtr<T> :: ListEmpty()
{
     while (head != NULL) {
           ListDelete(head);
     }
     length = 0;
}

template <class T> void ListPtr<T> :: DeleteData()
{
     Node<T*>* current = head;
     while (current != NULL)
     {
            delete current->data;
            current->data = NULL;
            current = current->next;
     }
     length = 0;
}

//-----------------------------------------------------------------------------

template <class T> Queue<T> :: Queue()
{
   ListQ = new List<T>();
}

template <class T> Queue<T> :: ~Queue()
{
   delete ListQ;
}

template <class T>
long int Queue<T> :: QueueLength()
{
   return ListQ->ListLength();
}

template <class T> bool Queue<T> :: Dequeue(T &rval)
{
   if (ListQ->ListLength() > 0) {
       rval = ListQ->head->data;
       ListQ->ListDelete(ListQ->head);
       return true;
   }
   else {
      return false;
   }
}

template <class T> void Queue<T> :: Enqueue(const T &item)
{
   ListQ->InsertRear(item);
}


template <class T> QueuePtr<T> :: QueuePtr()
{
   ListQ = new ListPtr<T>();
}

template <class T> QueuePtr<T> :: ~QueuePtr()
{
   delete ListQ;
}

template <class T>
long int QueuePtr<T> :: QueueLength()
{
   return ListQ->ListLength();
}

template <class T> bool QueuePtr<T> :: Dequeue(T * &rval)
{
   if (ListQ->ListLength() > 0) {
       rval = ListQ->head->data;
       ListQ->ListDelete(ListQ->head);
       return true;
   }
   else {
      return false;
   }
}

template <class T> void QueuePtr<T> :: Enqueue(T * &item)
{
   ListQ->InsertRear(item);
}

template <class T> void QueuePtr<T> :: DeleteData()
{
   ListQ->DeleteData();
}
	
		 
#endif

#ifndef STRASSEN_H
#define STRASSEN_H


struct elem_s {
  FLA_Obj *obj;
  FLA_Obj *objH;
  struct  elem_s *prev;
  struct  elem_s *next;
};

struct workspace_s {
  int size;
  struct elem_s *head;
  struct elem_s *tail;
};

typedef struct elem_s Strassen_Elem;
typedef struct workspace_s Strassen_Workspace;

Strassen_Workspace* Strassen_Workspace_new();
void Strassen_allocateObj(Strassen_Workspace *wks, FLA_Obj **objp, FLA_Obj **objHp, dim_t n, dim_t nb_alg);
void Strassen_releaseObj(Strassen_Elem *elem);
void Strassen_Workspace_free(Strassen_Workspace *wks);


FLA_Error FLA_Part_Even_2x2( FLA_Obj *A,  FLA_Obj *A11, FLA_Obj *A12,
                                    FLA_Obj *A21, FLA_Obj *A22);
 

FLA_Error FLA_Strassen( FLA_Obj *, FLA_Obj *, FLA_Obj *);

#endif

#include "FLAME.h"


struct elem_s {
  FLA_Obj *obj;
  FLA_Obj *objH;
  struct  elem_s *prev;
  struct  elem_s *next;
}

struct workspace_s {
  int size;
  struct obj_s *head;
  struct obj_s *tail;
}

typedef struct elem_s Strassen_Elem;
typedef struct workspace_s Strassen_Workspace;


void workspace_error(char* func_name, char* msg) {
  fprintf(stderr, "Workspace_Error: %s():  %s\n", func_name, msg);
}

Strassen_Workspace* Strassen_Workspace_new() {
  Strassen_Workspace *wks = (Strassen_Workspace *) malloc (sizeof(Strassen_Workspace));
  if (!wks) workspace_error("Workspace_new", "memory allocation failed");
  wks->head = NULL;
  wks->tail = NULL;
  wks->size = 0;
  return wks;
}

void Strassen_allocateObj(Strassen_Workspace *wks, FLA_Obj *obj, FLA_Obj *objH, int n) {
  obj = (FLA_Obj *)FLA_malloc(sizeof(FLA_Obj));
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, obj );
  FLA_Set( FLA_ZERO, obj );
  FLASH_Obj_create_without_buffer( FLA_Obj_datetype(obj), FLA_Obj_length(obj), FLA_Obj_width(obj), 1, &nb_alg, &objH );
  FLASH_attach_buffer(FLA_Obj_buffer_at_view(obj), FLA_Obj_row_stride(obj), FLA_Obj_col_stride(obj), &objH);
	
  Strassen_Elem *elem = (Strassen_Elem *)FLA_malloc(sizeof(Strassen_Elem));
  if(!elem) workspace_error("allocateObj", "memory allocation failed");
  elem->obj = obj;
  elem->objH = objH;
  elem->prev = NULL;
  elem->next = NULL;

  if (wks->size == 0) {
	wks->head = elem;
	wks->tail = elem;
  } else {
	wks->tail->next = elem;
	elem->prev = wks->tail;
	wks->tail = elem;
  }
  wks->size ++;
}

void Strassen_releaseObj(Strassen_Elem *elem) {
  FLA_Obj *obj = elem->obj, *objH = elem->objH;
  FLA_Obj_free( obj );
  FLA_free(obj);
  FLA_Obj_free( objH );
  FLA_free(objH);
  FLA_free(elem);
}


void Strassen_Workspace_free(Strassen_Workspace *wks) {
  Strassen_Elem *elem = wks->head, next;
  while (elem != NULL) {
	next = elem->next;
	Strassen_release(elem);
	elem = next;
  }
}







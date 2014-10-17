#include "FLAME.h"
#include "Strassen_prototypes.h"

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

void Strassen_allocateObj(Strassen_Workspace *wks, FLA_Obj *obj, FLA_Obj *objH, int n, int nb_alg) {
  Strassen_Elem *elem;
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, obj );
  FLA_Set( FLA_ZERO, *obj );
  FLASH_Obj_create_without_buffer( FLA_Obj_datatype(*obj), FLA_Obj_length(*obj), FLA_Obj_width(*obj), 1, &nb_alg, objH );
  //printf("n: %d\tnb_alg:%d\tobj_len:%d\tobj_width:%d\tobjH_len:%d\tobjH_wid:%d\n", n, nb_alg, FLA_Obj_length(*obj), FLA_Obj_width(*obj), FLA_Obj_length(*objH), FLA_Obj_width(*objH));
  FLASH_Obj_attach_buffer(FLA_Obj_buffer_at_view(*obj), FLA_Obj_row_stride(*obj), FLA_Obj_col_stride(*obj), objH);
	
  elem = (Strassen_Elem *)FLA_malloc(sizeof(Strassen_Elem));
  if(!elem) workspace_error("allocateObj", "memory allocation failed");
  elem->obj = obj;
  elem->objH = objH;
  elem->prev = NULL;
  elem->next = NULL;

  printf("elem->obj addr: %d\n", elem->obj);
  printf("obj: %d\n", obj);
  //printf("OBJdatatype:%d\n", FLA_Obj_datatype(*objH));
  //printf("wks->size:%d\n", wks->size);
  if (wks->size == 0) {
	printf("-------------------\n");
	wks->head = elem;
	wks->tail = elem;
	printf("head: %d\n", elem);
	printf("headelem->obj: %d\n", elem->obj);
	printf("headelem->objH: %d\n", elem->objH);
	printf("obj: %d\n", obj);
	printf("objH: %d\n", objH);
	//printf("elemtype:%d\n", FLA_Obj_elemtype(*(elem->obj)));
	//printf("m:%d\n", FLA_Obj_length(*(elem->obj)));
	//printf("elemtype:%d\n", FLA_Obj_elemtype(*(elem->objH)));

	FLA_Obj_show("obj:", *obj, "%11.3e", "...END...");

  } else {
	wks->tail->next = elem;
	elem->prev = wks->tail;
	wks->tail = elem;
  }
  wks->size ++;
}

void Strassen_releaseObj(Strassen_Elem *elem) {
  FLA_Obj *obj = elem->obj, *objH = elem->objH;

  //FLA_Obj_show("obj:", *obj, "%11.3e", "...END...");
  printf("elemtype:%d\n", FLA_Obj_elemtype(*(elem->obj)));
  printf("m:%d\n", FLA_Obj_length(*obj));
  printf("elemtype:%d\n", FLA_Obj_elemtype(*(elem->objH)));

  FLA_Obj_show("obj:", *obj, "%11.3e", "...END...");

  FLASH_Obj_free_without_buffer( objH );
  FLA_Obj_free( obj );
  FLA_free(elem);
}

void Strassen_Workspace_print(Strassen_Workspace *wks) {
  Strassen_Elem *elem = wks->head, *next;
  while (elem != NULL) {
	next = elem->next;

	printf("headelem->obj: %d\n", elem->obj);
	FLA_Obj_show("obj:", *(elem->obj), "%11.3e", "...END...");
	//FLA_Obj_show(*(elem->objH));
	elem = next;
  }
}


void Strassen_Workspace_free(Strassen_Workspace *wks) {
  Strassen_Elem *elem = wks->head, *next;
  printf("head: %d\n", elem);
  printf("headelem->obj: %d\n", elem->obj);
  printf("headelem->objH: %d\n", elem->objH);
  while (elem != NULL) {
	printf("elemtype:%d\n", FLA_Obj_elemtype(*(elem->obj)));
	next = elem->next;
	Strassen_releaseObj(elem);
	elem = next;
  }
}







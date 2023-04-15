/*
 * Copyright (c) 2023 Bryance Oyang
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *
 * === websocket_ctube ===
 *
 * View README.md for usage info and documentation.
 *
 * The contents of this file are copied from ./src by ./pkg.sh in the
 * websocket_ctube distribution.
 */

#ifndef WS_CTUBE_H
#define WS_CTUBE_H
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */



#ifndef WS_CTUBE_API_H
#define WS_CTUBE_API_H

#include <stddef.h>

#ifndef typeof
#define typeof __typeof__
#endif

struct ws_ctube;

/**
 * ws_ctube_open - create a ws_ctube websocket server. When finished, close with
 * ws_ctube_close()
 *
 * @param port port for websocket server
 * @param max_nclient maximum number of websocket client connections allowed
 * @param timeout_ms timeout (ms) for server starting and websocket handshake
 * or 0 for no timeout
 * @param max_broadcast_fps maximum number of broadcasts per second to rate
 * limit broadcasting or 0 for no limit. For best performance, disable by
 * setting to 0 and manually rate limit broadcasts.
 *
 * @return on success, a struct ws_ctube* is returned; on failure,
 * NULL is returned
 */
struct ws_ctube *ws_ctube_open(int port, int max_nclient, int timeout_ms, double max_broadcast_fps);

/**
 * ws_ctube_close - terminate ws_ctube server and cleanup
 */
void ws_ctube_close(struct ws_ctube *ctube);

/**
 * ws_ctube_broadcast - tries to queue data for sending to all connected
 * websocket clients.
 *
 * If max_broadcast_fps was nonzero when ws_ctube_open was called, this function
 * is rate-limited accordingly and returns failure if called too soon.
 *
 * Data is copied to an internal out-buffer, then this function returns. Actual
 * network operations will be handled internally and opaquely by separate
 * threads.
 *
 * Though non-blocking, try not to unnecessarily call this function in
 * performance-critical loops.
 *
 * If other threads can write to *data, get a read-lock to protect *data before
 * broadcasting. The read-lock can be released immediately once this function
 * returns.
 *
 * @param ctube the websocket ctube
 * @param data pointer to data to broadcast
 * @param data_size bytes of data
 *
 * @return 0 on success, nonzero otherwise
 */
int ws_ctube_broadcast(struct ws_ctube *ctube, const void *data, size_t data_size);

#endif /* WS_CTUBE_API_H */
#include <pthread.h>
#include <signal.h>
#include <stddef.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <time.h>
#include <stdint.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>


#ifndef WS_CTUBE_LIKELY_H
#define WS_CTUBE_LIKELY_H

#define ws_ctube_likely(x) __builtin_expect(!!(x), 1)
#define ws_ctube_unlikely(x) __builtin_expect(!!(x), 0)

#endif /* WS_CTUBE_LIKELY_H */
#ifndef WS_CTUBE_CONTAINER_OF_H
#define WS_CTUBE_CONTAINER_OF_H

/**
 * get the container of ptr which is a member of a struct
 *
 * Example:
 *
 * struct a {
 * 	struct b bmember;
 * 	struct c cmember;
 * };
 *
 * struct a *a_ptr;
 * struct c *c_ptr = &a_ptr->cmember;
 *
 * assert(a_ptr == ws_ctube_container_of(c_ptr, struct a, cmember))
 */
#define ws_ctube_container_of(ptr, type, member) ({ \
	typeof(((type *)0)->member) *_ws_ctube_container_of_ptr = (ptr); \
	((type *)((char *)(_ws_ctube_container_of_ptr) - offsetof(type, member)));})

#endif /* WS_CTUBE_CONTAINER_OF_H */




#ifndef WS_CTUBE_REF_COUNT_H
#define WS_CTUBE_REF_COUNT_H


/** including this in a larger struct allows it to be reference counted */
struct ws_ctube_ref_count {
	volatile int refc;
};

static int ws_ctube_ref_count_init(struct ws_ctube_ref_count *ref_count)
{
	__atomic_store_n(&ref_count->refc, (int)0, __ATOMIC_SEQ_CST);
	return 0;
}

static void ws_ctube_ref_count_destroy(struct ws_ctube_ref_count *ref_count)
{
	__atomic_store_n(&ref_count->refc, (int)(-1), __ATOMIC_SEQ_CST);
}

#ifndef __cplusplus
#define ws_ctube_ref_count_acquire(ptr, ref_count_member) do { \
		_Static_assert(__builtin_types_compatible_p(typeof((ptr)->ref_count_member), struct ws_ctube_ref_count), "type mismatch in ws_ctube_ref_count_acquire()"); \
		__atomic_add_fetch(&(ptr)->ref_count_member.refc, (int)1, __ATOMIC_SEQ_CST); \
	} while (0);

#define ws_ctube_ref_count_release(ptr, ref_count_member, release_routine) do { \
		_Static_assert(__builtin_types_compatible_p(typeof((ptr)->ref_count_member), struct ws_ctube_ref_count), "type mismatch in ws_ctube_ref_count_release()"); \
		const int _ws_ctube_ref_count_refc = __atomic_sub_fetch(&(ptr)->ref_count_member.refc, (int)1, __ATOMIC_SEQ_CST); \
		if (_ws_ctube_ref_count_refc <= 0) { \
			if (ws_ctube_likely(_ws_ctube_ref_count_refc == 0)) { \
				release_routine(ptr); \
			} else { \
				raise(SIGSEGV); \
			} \
		} \
	} while (0);
#else /* __cplusplus */
#define ws_ctube_ref_count_acquire(ptr, ref_count_member) do { \
		__atomic_add_fetch(&(ptr)->ref_count_member.refc, (int)1, __ATOMIC_SEQ_CST); \
	} while (0);

#define ws_ctube_ref_count_release(ptr, ref_count_member, release_routine) do { \
		const int _ws_ctube_ref_count_refc = __atomic_sub_fetch(&(ptr)->ref_count_member.refc, (int)1, __ATOMIC_SEQ_CST); \
		if (_ws_ctube_ref_count_refc <= 0) { \
			if (ws_ctube_likely(_ws_ctube_ref_count_refc == 0)) { \
				release_routine(ptr); \
			} else { \
				raise(SIGSEGV); \
			} \
		} \
	} while (0);
#endif /* __cplusplus */

#endif /* WS_CTUBE_REF_COUNT_H */




#ifndef WS_CTUBE_LIST_H
#define WS_CTUBE_LIST_H


/** including this in a larger struct allows it to be a part of a list */
struct ws_ctube_list_node {
	struct ws_ctube_list_node *prev;
	struct ws_ctube_list_node *next;
	pthread_mutex_t mutex;
};

static int ws_ctube_list_node_init(struct ws_ctube_list_node *node)
{
	node->prev = NULL;
	node->next = NULL;
	pthread_mutex_init(&node->mutex, NULL);
	return 0;
}

static void ws_ctube_list_node_destroy(struct ws_ctube_list_node *node)
{
	node->prev = NULL;
	node->next = NULL;
	pthread_mutex_destroy(&node->mutex);
}

/** thread-safe circular doubly-linked list */
struct ws_ctube_list {
	struct ws_ctube_list_node head;
	int len;
	pthread_mutex_t mutex;
};

static int ws_ctube_list_init(struct ws_ctube_list *l)
{
	ws_ctube_list_node_init(&l->head);
	l->head.next = &l->head;
	l->head.prev = &l->head;
	l->len = 0;
	pthread_mutex_init(&l->mutex, NULL);
	return 0;
}

static void ws_ctube_list_destroy(struct ws_ctube_list *l)
{
	ws_ctube_list_node_destroy(&l->head);
	pthread_mutex_destroy(&l->mutex);
}

static void _ws_ctube_list_add_after(struct ws_ctube_list_node *prior, struct ws_ctube_list_node *new_node)
{
	new_node->prev = prior;
	new_node->next = prior->next;
	prior->next->prev = new_node;
	prior->next = new_node;
}

static void _ws_ctube_list_node_unlink(struct ws_ctube_list_node *node)
{
	node->next->prev = node->prev;
	node->prev->next = node->next;
	node->next = NULL;
	node->prev = NULL;
}

static inline void ws_ctube_list_unlink(struct ws_ctube_list *l, struct ws_ctube_list_node *node)
{
	pthread_mutex_lock(&l->mutex);
	pthread_mutex_lock(&node->mutex);
	_ws_ctube_list_node_unlink(node);
	l->len--;
	pthread_mutex_unlock(&node->mutex);
	pthread_mutex_unlock(&l->mutex);
}

static inline int ws_ctube_list_push_front(struct ws_ctube_list *l, struct ws_ctube_list_node *node)
{
	int retval = 0;
	pthread_mutex_lock(&l->mutex);
	pthread_mutex_lock(&node->mutex);

	if (node->next != NULL || node->prev != NULL) {
		retval = -1;
		goto out;
	}
	_ws_ctube_list_add_after(&l->head, node);
	l->len++;

out:
	pthread_mutex_unlock(&node->mutex);
	pthread_mutex_unlock(&l->mutex);
	return retval;
}

static inline int ws_ctube_list_push_back(struct ws_ctube_list *l, struct ws_ctube_list_node *node)
{
	int retval = 0;
	pthread_mutex_lock(&l->mutex);
	pthread_mutex_lock(&node->mutex);

	if (node->next != NULL || node->prev != NULL) {
		retval = -1;
		goto out;
	}
	_ws_ctube_list_add_after(l->head.prev, node);
	l->len++;

out:
	pthread_mutex_unlock(&node->mutex);
	pthread_mutex_unlock(&l->mutex);
	return retval;
}

static inline struct ws_ctube_list_node *ws_ctube_list_pop_front(struct ws_ctube_list *l)
{
	struct ws_ctube_list_node *front;

	pthread_mutex_lock(&l->mutex);
	if (l->len == 0) {
		pthread_mutex_unlock(&l->mutex);
		return NULL;
	}
	front = l->head.next;
	pthread_mutex_lock(&front->mutex);

	_ws_ctube_list_node_unlink(front);
	l->len--;

	pthread_mutex_unlock(&front->mutex);
	pthread_mutex_unlock(&l->mutex);
	return front;
}

static inline struct ws_ctube_list_node *ws_ctube_list_pop_back(struct ws_ctube_list *l)
{
	struct ws_ctube_list_node *back;

	pthread_mutex_lock(&l->mutex);
	if (l->len == 0) {
		pthread_mutex_unlock(&l->mutex);
		return NULL;
	}
	back = l->head.prev;
	pthread_mutex_lock(&back->mutex);

	_ws_ctube_list_node_unlink(back);
	l->len--;

	pthread_mutex_unlock(&back->mutex);
	pthread_mutex_unlock(&l->mutex);
	return back;
}

/** loop through list_node */
#define ws_ctube_list_for_each(list, node) \
	for ( \
	node = (list)->head.next; \
	node != &((list)->head); \
	node = node->next)

/** loop through containers of list_node */
#define ws_ctube_list_for_each_entry(list, entry, member) \
	for ( \
	entry = (list)->head.next != &((list)->head) ? \
	ws_ctube_container_of((list)->head.next, typeof(*entry), member) : \
	NULL; \
	entry != NULL; \
	entry = entry->member.next != &((list)->head) ? \
	ws_ctube_container_of(entry->member.next, typeof(*entry), member) : \
	NULL)

#endif /* WS_CTUBE_LIST_H */


#ifndef WS_CTUBE_CRYPT_H
#define WS_CTUBE_CRYPT_H


void ws_ctube_b64_encode(unsigned char *out, const unsigned char *in, size_t in_bytes);
void ws_ctube_sha1sum(unsigned char *out, const unsigned char *in, size_t len);

#endif /* WS_CTUBE_CRYPT_H */




#ifndef WS_CTUBE_SOCKET_H
#define WS_CTUBE_SOCKET_H


#ifndef MSG_NOSIGNAL
#define MSG_NOSIGNAL 0
#endif

/**
 * send all characters in buf through a socket
 *
 * @param fd file descriptor
 * @param buf buffer
 * @param buf_size size of buffer in bytes
 *
 * @return 0 on success, -1 otherwise
 */
static inline int ws_ctube_socket_send_all(const int fd, const char *buf, ssize_t buf_size)
{
	while (buf_size > 0) {
		ssize_t nsent = send(fd, buf, buf_size, MSG_NOSIGNAL);
		if (nsent < 1) {
			return -1;
		}
		buf += nsent;
		buf_size -= nsent;
	}

	return 0;
}

/**
 * receive all characters up to buf_size or when delim is encountered
 *
 * @param fd file descriptor
 * @param buf buffer to read into
 * @param buf_size size of buffer (including '\0' if delim is not NULL)
 * @param delim should be a null-terminated string or NULL
 * if delim is NULL, receive up to buf_size bytes
 * if delim is not NULL, receive up to buf_size - 1 bytes or until delim is found
 *
 * @return 0 on success, -1 otherwise
 */
static inline int ws_ctube_socket_recv_all(const int fd, char *const buf, const ssize_t buf_size, const char *delim)
{
	char *remain_buf = buf;
	ssize_t remain_size;

	if (buf_size <= 0) {
		return 0;
	}

	if (delim != NULL) {
		remain_size = buf_size - 1;
		buf[buf_size - 1] = '\0';
	} else {
		remain_size = buf_size;
	}

	while (remain_size > 0) {
		ssize_t nrecv = recv(fd, remain_buf, remain_size, MSG_NOSIGNAL);
		if (nrecv < 1) {
			return -1;
		}
		remain_buf += nrecv;
		remain_size -= nrecv;

		/* look for delim */
		if (delim != NULL && remain_size > 0) {
			*remain_buf = '\0';
			char *dpos = strstr(buf, delim);
			if (dpos == NULL) {
				continue;
			}

			dpos += strlen(delim);
			*dpos = '\0';
			return 0;
		}
	}

	return 0;
}

static inline int ws_ctube_bind_server(int server_sock, int port)
{
	struct sockaddr_in sa;
	sa.sin_family = AF_INET;
	sa.sin_addr.s_addr = htonl(INADDR_ANY);
	sa.sin_port = htons(port);
	return bind(server_sock, (struct sockaddr *)&sa, sizeof(sa));
}

#endif /* WS_CTUBE_SOCKET_H */


#ifndef WS_CTUBE_WS_BASE_H
#define WS_CTUBE_WS_BASE_H


#define WS_CTUBE_FRAME_HDR_SIZE 2
#define WS_CTUBE_MAX_PAYLD_SIZE 125

/**
 * make a websocket frame
 *
 * @param frame pointer to buffer where frame shall be written; needs to have
 * size of at least WS_CTUBE_FRAME_HDR_SIZE + WS_CTUBE_MAX_PAYLD_SIZE bytes
 * @param msg pointer to data
 * @param msg_size bytes of message
 * @param first whether this is the first frame in a sequence
 *
 * @return number of bytes of msg contained in frame
 */
int ws_ctube_ws_mkframe(char *frame, const char *msg, size_t msg_size, int first);

int ws_ctube_ws_send(int conn, const char *msg, size_t msg_size);
int ws_ctube_ws_recv(int conn, char *msg, int *msg_size, size_t max_msg_size);
int ws_ctube_ws_is_ping(const char *msg, int msg_size);
int ws_ctube_ws_pong(int conn, const char *msg, int msg_size);
int ws_ctube_ws_handshake(int conn, const struct timeval *timeout);

#endif /* WS_CTUBE_WS_BASE_H */


#ifndef WS_CTUBE_STRUCT_H
#define WS_CTUBE_STRUCT_H


/** holds data to be sent/received over the network */
struct ws_ctube_data {
	void *data;
	size_t data_size;

	pthread_mutex_t mutex;
	struct ws_ctube_list_node lnode;
	struct ws_ctube_ref_count refc;
};

static int ws_ctube_data_init(struct ws_ctube_data *ws_ctube_data, const void *data, size_t data_size)
{
	if (data_size > 0) {
		ws_ctube_data->data = (typeof(ws_ctube_data->data))malloc(data_size);
		if (ws_ctube_data->data == NULL) {
			goto out_nodata;
		}

		if (data != NULL) {
			memcpy(ws_ctube_data->data, data, data_size);
		}
	}

	ws_ctube_data->data_size = data_size;

	pthread_mutex_init(&ws_ctube_data->mutex, NULL);
	ws_ctube_list_node_init(&ws_ctube_data->lnode);
	ws_ctube_ref_count_init(&ws_ctube_data->refc);
	return 0;

out_nodata:
	return -1;
}

static void ws_ctube_data_destroy(struct ws_ctube_data *ws_ctube_data)
{
	if (ws_ctube_data->data != NULL) {
		free(ws_ctube_data->data);
		ws_ctube_data->data = NULL;
	}

	ws_ctube_data->data_size = 0;

	pthread_mutex_destroy(&ws_ctube_data->mutex);
	ws_ctube_list_node_destroy(&ws_ctube_data->lnode);
	ws_ctube_ref_count_destroy(&ws_ctube_data->refc);
}

/** copy data into ws_ctube_data */
static inline int ws_ctube_data_cp(struct ws_ctube_data *ws_ctube_data, const void *data, size_t data_size)
{
	int retval = 0;

	pthread_mutex_lock(&ws_ctube_data->mutex);
	if (ws_ctube_data->data_size < data_size) {
		if (ws_ctube_data->data != NULL) {
			free(ws_ctube_data->data);
		}

		ws_ctube_data->data = (typeof(ws_ctube_data->data))malloc(data_size);
		if (ws_ctube_data->data == NULL) {
			retval = -1;
			goto out;
		}

		ws_ctube_data->data_size = data_size;
	}

	memcpy(ws_ctube_data->data, data, data_size);

out:
	pthread_mutex_unlock(&ws_ctube_data->mutex);
	return retval;
}

static void ws_ctube_data_free(struct ws_ctube_data *ws_ctube_data)
{
	ws_ctube_data_destroy(ws_ctube_data);
	free(ws_ctube_data);
}

/** represents a client connection and owns their associated reader/writer threads */
struct ws_ctube_conn_struct {
	int fd;
	struct ws_ctube *ctube;

	/* to prevent double shutdown */
	int stopping;
	pthread_mutex_t stopping_mutex;

	/** reader thread */
	pthread_t reader_tid;
	/** writer thread */
	pthread_t writer_tid;

	struct ws_ctube_ref_count refc;
	struct ws_ctube_list_node lnode;
};

static int ws_ctube_conn_struct_init(struct ws_ctube_conn_struct *conn, int fd, struct ws_ctube *ctube)
{
	conn->fd = fd;
	conn->ctube = ctube;

	conn->stopping = 0;
	pthread_mutex_init(&conn->stopping_mutex, NULL);

	ws_ctube_ref_count_init(&conn->refc);
	ws_ctube_list_node_init(&conn->lnode);
	return 0;
}

static void ws_ctube_conn_struct_destroy(struct ws_ctube_conn_struct *conn)
{
	int oldstate, statevar;
	pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &oldstate);
	close(conn->fd);
	pthread_setcancelstate(oldstate, &statevar);

	conn->fd = -1;
	conn->ctube = NULL;

	conn->stopping = 0;
	pthread_mutex_destroy(&conn->stopping_mutex);

	ws_ctube_ref_count_destroy(&conn->refc);
	ws_ctube_list_node_destroy(&conn->lnode);
}

static void ws_ctube_conn_struct_free(struct ws_ctube_conn_struct *conn)
{
	ws_ctube_conn_struct_destroy(conn);
	free(conn);
}

enum ws_ctube_qaction {
	WS_CTUBE_CONN_START,
	WS_CTUBE_CONN_STOP
};

/** a work item in FIFO queue: start or stop a connection to client */
struct ws_ctube_conn_qentry {
	struct ws_ctube_conn_struct *conn;
	enum ws_ctube_qaction act;
	struct ws_ctube_list_node lnode;
};

static int ws_ctube_conn_qentry_init(struct ws_ctube_conn_qentry *qentry, struct ws_ctube_conn_struct *conn, enum ws_ctube_qaction act)
{
	ws_ctube_ref_count_acquire(conn, refc);
	qentry->conn = conn;
	qentry->act = act;
	ws_ctube_list_node_init(&qentry->lnode);
	return 0;
}

static void ws_ctube_conn_qentry_destroy(struct ws_ctube_conn_qentry *qentry)
{
	ws_ctube_ref_count_release(qentry->conn, refc, ws_ctube_conn_struct_free);
	qentry->conn = NULL;
	ws_ctube_list_node_destroy(&qentry->lnode);
}

static void ws_ctube_conn_qentry_free(struct ws_ctube_conn_qentry *qentry)
{
	ws_ctube_conn_qentry_destroy(qentry);
	free(qentry);
}

/** main struct for ws_ctube */
struct ws_ctube {
	int server_sock;
	int port;
	int max_nclient;

	/* to have timeout on operations */
	struct timespec timeout_spec;
	struct timeval timeout_val;

	/* currently unused (for incoming data) */
	struct ws_ctube_list in_data_list;
	int in_data_pred;
	pthread_mutex_t in_data_mutex;
	pthread_cond_t in_data_cond;

	/* current ws_ctube_data representing data to be sent */
	struct ws_ctube_data *out_data;
	unsigned long out_data_id;
	pthread_mutex_t out_data_mutex;
	pthread_cond_t out_data_cond;

	/* rate-limit broadcasting */
	double max_bcast_fps;
	struct timespec prev_bcast_time;

	/* the FIFO work queue: connection handler starts/stops client
	 * connections based on queued actions */
	struct ws_ctube_list connq;
	int connq_pred;
	pthread_mutex_t connq_mutex;
	pthread_cond_t connq_cond;

	/* allows ws_ctube_open() to know if server successfully started or not */
	int server_inited;
	pthread_mutex_t server_init_mutex;
	pthread_cond_t server_init_cond;

	/** connection handler thread */
	pthread_t handler_tid;
	/** server thread */
	pthread_t server_tid;
};

static int ws_ctube_init(
	struct ws_ctube *ctube,
	int port,
	int max_nclient,
	unsigned int timeout_ms,
	double max_broadcast_fps)
{
	ctube->server_sock = -1;
	ctube->port = port;
	ctube->max_nclient = max_nclient;

	ctube->timeout_spec.tv_sec = timeout_ms / 1000;
	ctube->timeout_spec.tv_nsec = (timeout_ms % 1000) * 1000000;
	ctube->timeout_val.tv_sec = timeout_ms / 1000;
	ctube->timeout_val.tv_usec = (timeout_ms % 1000) * 1000;

	ws_ctube_list_init(&ctube->in_data_list);
	ctube->in_data_pred = 0;
	pthread_mutex_init(&ctube->in_data_mutex, NULL);
	pthread_cond_init(&ctube->in_data_cond, NULL);

	ctube->out_data = NULL;
	ctube->out_data_id = 0;
	pthread_mutex_init(&ctube->out_data_mutex, NULL);
	pthread_cond_init(&ctube->out_data_cond, NULL);

	ctube->max_bcast_fps = max_broadcast_fps;
	ctube->prev_bcast_time.tv_sec = 0;
	ctube->prev_bcast_time.tv_nsec = 0;

	ws_ctube_list_init(&ctube->connq);
	ctube->connq_pred = 0;
	pthread_mutex_init(&ctube->connq_mutex, NULL);
	pthread_cond_init(&ctube->connq_cond, NULL);

	ctube->server_inited = 0;
	pthread_mutex_init(&ctube->server_init_mutex, NULL);
	pthread_cond_init(&ctube->server_init_cond, NULL);

	return 0;
}

static void _ws_ctube_data_list_clear(struct ws_ctube_list *dlist)
{
	struct ws_ctube_list_node *node;
	struct ws_ctube_data *data;

	while ((node = ws_ctube_list_pop_front(dlist)) != NULL) {
		data = ws_ctube_container_of(node, typeof(*data), lnode);
		ws_ctube_data_free(data);
	}
}

static void _ws_ctube_connq_clear(struct ws_ctube_list *connq)
{
	struct ws_ctube_list_node *node;
	struct ws_ctube_conn_qentry *qentry;

	while ((node = ws_ctube_list_pop_front(connq)) != NULL) {
		qentry = ws_ctube_container_of(node, typeof(*qentry), lnode);
		ws_ctube_conn_qentry_free(qentry);
	}
}

static void ws_ctube_destroy(struct ws_ctube *ctube)
{
	ctube->server_sock = -1;
	ctube->port = -1;
	ctube->max_nclient = -1;

	ctube->timeout_spec.tv_sec = 0;
	ctube->timeout_spec.tv_nsec = 0;
	ctube->timeout_val.tv_sec = 0;
	ctube->timeout_val.tv_usec = 0;

	_ws_ctube_data_list_clear(&ctube->in_data_list);
	ws_ctube_list_destroy(&ctube->in_data_list);
	ctube->in_data_pred = 0;
	pthread_mutex_destroy(&ctube->in_data_mutex);
	pthread_cond_destroy(&ctube->in_data_cond);

	if (ctube->out_data != NULL) {
		ws_ctube_ref_count_release(ctube->out_data, refc, ws_ctube_data_free);
		ctube->out_data = NULL;
	}
	ctube->out_data_id = 0;
	pthread_mutex_destroy(&ctube->out_data_mutex);
	pthread_cond_destroy(&ctube->out_data_cond);

	ctube->max_bcast_fps = 0;
	ctube->prev_bcast_time.tv_sec = 0;
	ctube->prev_bcast_time.tv_nsec = 0;

	_ws_ctube_connq_clear(&ctube->connq);
	ws_ctube_list_destroy(&ctube->connq);
	ctube->connq_pred = 0;
	pthread_mutex_destroy(&ctube->connq_mutex);
	pthread_cond_destroy(&ctube->connq_cond);

	ctube->server_inited = 0;
	pthread_mutex_destroy(&ctube->server_init_mutex);
	pthread_cond_destroy(&ctube->server_init_cond);
}

#endif /* WS_CTUBE_STRUCT_H */






static volatile int ws_ctube_b64_encode_inited = 0;
static volatile unsigned char ws_ctube_b64_encode_table[64];
static const uint32_t ws_ctube_b64_mask = 63;

/**
 * 0-25: A-Z
 * 26-51: a-z
 * 52-61: 0-9
 * 62: +
 * 63: /
 */
static void ws_ctube_init_b64_encode_table()
{
	/* prevent double init */
	if (__atomic_exchange_n(&ws_ctube_b64_encode_inited, (int)1, __ATOMIC_SEQ_CST))
		return;

	for (int i = 0; i < 26; i++) {
		ws_ctube_b64_encode_table[i] = 'A' + i;
		ws_ctube_b64_encode_table[26 + i] = 'a' + i;
	}
	for (int i = 0; i < 10; i++) {
		ws_ctube_b64_encode_table[52 + i] = '0' + i;
	}
	ws_ctube_b64_encode_table[62] = '+';
	ws_ctube_b64_encode_table[63] = '/';
}

/* 3*8 -> 4*6 base64 encoding */
static void ws_ctube_b64_encode_triplet(unsigned char *out, unsigned char in0, unsigned char in1, unsigned char in2)
{
	uint32_t triplet = ((uint32_t)in0 << 16) + ((uint32_t)in1 << 8) + (uint32_t)in2;
	for (int i = 0; i < 4; i++) {
		uint32_t b64_val = (triplet >> (6*(3 - i))) & ws_ctube_b64_mask;
		out[i] = ws_ctube_b64_encode_table[b64_val];
	}
}

/**
 * gives the base64 encoding of in. out must be large enough to hold output + 1
 * for null terminator
 */
void ws_ctube_b64_encode(unsigned char *out, const unsigned char *in, size_t in_bytes)
{
#ifndef __cplusplus
	_Static_assert(sizeof(uint8_t) == sizeof(unsigned char), "ws_ctube_b64_encode(): unsigned char not 8 bits");
#endif /* __cplusplus */

	ws_ctube_init_b64_encode_table();

	while (in_bytes >= 3) {
		ws_ctube_b64_encode_triplet(out, in[0], in[1], in[2]);
		in_bytes -= 3;
		in += 3;
		out += 4;
	}

	if (in_bytes == 0) {
		out[0] = '\0';
	} else if (in_bytes == 1) {
		ws_ctube_b64_encode_triplet(out, in[0], 0, 0);
		out[2] = '=';
		out[3] = '=';
		out[4] = '\0';
	} else if (in_bytes == 2) {
		ws_ctube_b64_encode_triplet(out, in[0], in[1], 0);
		out[3] = '=';
		out[4] = '\0';
	}
}

static inline uint32_t ws_ctube_left_rotate(uint32_t w, uint32_t amount)
{
	return (w << amount) | (w >> (32 - amount));
}

/** see sha1 standard */
static inline void ws_ctube_sha1_wordgen(uint32_t *const word)
{
	for (int i = 16; i < 80; i++) {
		word[i] = ws_ctube_left_rotate(word[i-3] ^ word[i-8] ^ word[i-14] ^ word[i-16], 1);
	}
}

/** generate words 14,15 from input length (see sha1 standard) */
static inline void ws_ctube_sha1_total_len_cpy(uint32_t *words, const uint64_t total_bits)
{
	const uint64_t mask = 0xFF;

	words[14] = 0;
	words[15] = 0;
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 4; j++) {
			uint32_t byte = (total_bits >> 8*(8 - (4*i+j) - 1)) & mask;
			words[14 + i] |= byte << 8*(4 - j - 1);
		}
	}
}

/** copy byte stream into sha1 32-byte words */
static inline void ws_ctube_sha1_cp_to_words(uint32_t *const words, const uint8_t *const in, const size_t len, int pad)
{
	size_t i, w_ind, byte_ind;

	for (i = 0, w_ind = 0, byte_ind = 0; i < len; i++, w_ind = i / 4, byte_ind = i % 4) {
		uint32_t tmp = in[i];
		words[w_ind] |= tmp << 8*(4 - byte_ind - 1);
	}
	if (pad) {
		uint32_t one = 0x80;
		words[w_ind] |= one << 8*(4 - byte_ind - 1);
	}
}

/** genearte all sha1 words (see sha1 standard) */
static inline void ws_ctube_sha1_mkwords(uint32_t *const words, const uint8_t *const in, const size_t len, const int mode, const uint64_t total_bits)
{
	for (int i = 0; i < 16; i++) {
		words[i] = 0;
	}

	switch (mode) {
	case 0:
		if (len < 56) {
			ws_ctube_sha1_cp_to_words(words, in, len, 1);
			ws_ctube_sha1_total_len_cpy(words, total_bits);
		} else if (len < 64) {
			ws_ctube_sha1_cp_to_words(words, in, len, 1);
		} else {
			ws_ctube_sha1_cp_to_words(words, in, 64, 0);
		}
		break;

	case 1:
		/* needs total_len appended only */
		ws_ctube_sha1_total_len_cpy(words, total_bits);
		break;
	}

	ws_ctube_sha1_wordgen(words);
}

/**
 * computes the sha1 hash of in and stores the result in out. out must be able
 * to hold the 20 byte output of sha1
 */
void ws_ctube_sha1sum(unsigned char *out, const unsigned char *in, size_t len_bytes)
{
#ifndef __cplusplus
	_Static_assert(sizeof(uint8_t) == sizeof(unsigned char), "ws_ctube_sha1sum(): unsigned char not 8 bits");
#endif /* __cplusplus */

	const uint8_t *in_byte = (uint8_t *)in;
	uint8_t *const out_byte = (uint8_t *)out;
	const uint64_t total_bits = ((uint64_t)8 * (uint64_t)len_bytes);

	const uint32_t K[4] = {
		0x5A827999,
		0x6ED9EBA1,
		0x8F1BBCDC,
		0xCA62C1D6
	};

	uint32_t h[5] = {
		0x67452301,
		0xEFCDAB89,
		0x98BADCFE,
		0x10325476,
		0xC3D2E1F0
	};

	uint32_t words[80];
	int mode = 0;
	while(1) {
		ws_ctube_sha1_mkwords(words, in_byte, len_bytes, mode, total_bits);

		uint32_t A = h[0];
		uint32_t B = h[1];
		uint32_t C = h[2];
		uint32_t D = h[3];
		uint32_t E = h[4];

		for (int i = 0; i < 20; i++) {
			uint32_t temp = ws_ctube_left_rotate(A, 5) + ((B & C) | ((~B) & D)) + E + words[i] + K[0];
			E = D;
			D = C;
			C = ws_ctube_left_rotate(B, 30);
			B = A;
			A = temp;
		}
		for (int i = 20; i < 40; i++) {
			uint32_t temp = ws_ctube_left_rotate(A, 5) + (B ^ C ^ D) + E + words[i] + K[1];
			E = D;
			D = C;
			C = ws_ctube_left_rotate(B, 30);
			B = A;
			A = temp;
		}
		for (int i = 40; i < 60; i++) {
			uint32_t temp = ws_ctube_left_rotate(A, 5) + ((B & C) | (B & D) | (C & D)) + E + words[i] + K[2];
			E = D;
			D = C;
			C = ws_ctube_left_rotate(B, 30);
			B = A;
			A = temp;
		}
		for (int i = 60; i < 80; i++) {
			uint32_t temp = ws_ctube_left_rotate(A, 5) + (B ^ C ^ D) + E + words[i] + K[3];
			E = D;
			D = C;
			C = ws_ctube_left_rotate(B, 30);
			B = A;
			A = temp;
		}

		h[0] += A;
		h[1] += B;
		h[2] += C;
		h[3] += D;
		h[4] += E;

		if (mode != 0) {
			break;
		}
		if (len_bytes < 56) {
			break;
		} else if (len_bytes < 64) {
			mode = 1;
		} else {
			len_bytes -= 64;
			in_byte += 64;
		}
	}

	uint32_t mask = 0xFF;
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 4; j++) {
			out_byte[4*i + j] = (h[i] >> 8*(4 - j - 1)) & mask;
		}
	}
}







#define WS_DEBUG 0
#define WS_BUFLEN 4096

static void ws_print_frame(const char *prefix, char *frame, int len)
{
	if (!WS_DEBUG) {
		return;
	}

	printf("%s\n", prefix);
	for (int i = 0; i < len; i++) {
		for (int j = 7; j >= 0; j--) {
			uint8_t mask = ((uint8_t)1 << j);
			uint8_t val = ((frame[i] & mask) >> j);
			printf("%u", val);
		}
		printf("|");
		if (i % 4 == 3) {
			printf("\n");
		}
	}
	printf("\n\n");
	fflush(stdout);
}

/** create data frame according to websocket standard */
int ws_ctube_ws_mkframe(char *frame, const char *msg, size_t msg_size, int first)
{
	int payld_size;

	first = !!first;
	if (msg_size > WS_CTUBE_MAX_PAYLD_SIZE) {
		frame[0] = 2*first;
		payld_size = WS_CTUBE_MAX_PAYLD_SIZE;
	} else {
		frame[0] = 0b10000000 + 2*first;
		payld_size = msg_size;
	}

	frame[1] = payld_size;
	memcpy(&frame[WS_CTUBE_FRAME_HDR_SIZE], msg, payld_size);

	return payld_size;
}

/** send data in data frames according to websocket standard */
int ws_ctube_ws_send(int conn, const char *msg, size_t msg_size)
{
	int payld_size;
	char frame[WS_BUFLEN];

	for (int first = 1; msg_size > 0; first = 0, msg += payld_size, msg_size -= payld_size) {
		payld_size = ws_ctube_ws_mkframe(frame, msg, msg_size, first);
		const int frame_len = payld_size + WS_CTUBE_FRAME_HDR_SIZE;
		ws_print_frame("ws_ctube_ws_send()", frame, frame_len);
		if (ws_ctube_socket_send_all(conn, frame, frame_len) != 0) {
			return -1;
		}
	}

	return 0;
}

int ws_ctube_ws_recv(int conn, char *msg, int *msg_size, size_t max_msg_size)
{
	/* TODO */
	(void)conn;
	(void)msg;
	(void)msg_size;
	(void)max_msg_size;
	return -1;
}

int ws_ctube_ws_is_ping(const char *msg, int msg_size)
{
	/* TODO */
	(void)msg;
	(void)msg_size;
	return 0;
}

int ws_ctube_ws_pong(int conn, const char *msg, int msg_size)
{
	/* TODO */
	(void)conn;
	(void)msg;
	(void)msg_size;
	return 0;
}

/** extract the client key from handshake */
static char *ws_client_key(char *rbuf)
{
	char *client_key, *client_key_end;

	client_key = strstr(rbuf, "Sec-WebSocket-Key: ");
	if (client_key == NULL) {
		return NULL;
	}

	client_key += strlen("Sec-WebSocket-Key: ");
	client_key_end = strstr(client_key, "\r");
	if (client_key_end == NULL) {
		return NULL;
	}

	*client_key_end = '\0';
	if (WS_DEBUG) {
		printf("wskey\n%s\n", client_key);
	}
	return client_key;
}

/** compute server response key per websocket standard */
static int ws_server_response_key(char *server_key, const char *client_key)
{
	size_t sha1_hash_len = 20;
	char magic_client_key[WS_BUFLEN];
	char client_sha1[WS_BUFLEN];

	if (snprintf(magic_client_key, WS_BUFLEN, "%s258EAFA5-E914-47DA-95CA-C5AB0DC85B11", client_key) >= WS_BUFLEN) {
		return -1;
	}

	ws_ctube_sha1sum((unsigned char *)client_sha1, (unsigned char *)magic_client_key, strlen(magic_client_key));
	ws_ctube_b64_encode((unsigned char *)server_key, (unsigned char *)client_sha1, sha1_hash_len);

	return 0;
}

int ws_ctube_ws_handshake(int conn, const struct timeval *timeout)
{
	char rbuf[WS_BUFLEN];
	char *client_key;
	char server_key[WS_BUFLEN];
	char response[2*WS_BUFLEN];

	const char *const response_fmt = "HTTP/1.1 101 Switching Protocols\r\n"
				"Upgrade: websocket\r\n"
				"Connection: Upgrade\r\n"
				"Sec-WebSocket-Accept: %s\r\n\r\n";

	/* receive with timeout, but reset to old timeout afterwards */
	struct timeval old_timeout;
	socklen_t timeval_size = sizeof(old_timeout);
	if (getsockopt(conn, SOL_SOCKET, SO_RCVTIMEO, &old_timeout, &timeval_size) < 0) {
		goto err;
	}
	if (setsockopt(conn, SOL_SOCKET, SO_RCVTIMEO, timeout, sizeof(*timeout)) < 0) {
		goto err;
	}
	if (ws_ctube_socket_recv_all(conn, rbuf, WS_BUFLEN, "\r\n\r\n") != 0) {
		goto err;
	}
	if (setsockopt(conn, SOL_SOCKET, SO_RCVTIMEO, &old_timeout, sizeof(old_timeout)) < 0) {
		goto err;
	}

	/* ensure null termination of received data */
	rbuf[WS_BUFLEN - 1] = '\0';

	if (WS_DEBUG) {
		printf("get\n%s\n", rbuf);
	}

	client_key = ws_client_key(rbuf);
	if (client_key == NULL) {
		goto err;
	}
	if (ws_server_response_key(server_key, client_key) != 0) {
		goto err;
	}

	snprintf(response, sizeof(response)/sizeof(response[0]), response_fmt, server_key);
	if (WS_DEBUG) {
		printf("server response\n%s\n", response);
	}

	/* send with timeout, but reset to old timeout afterwards */
	if (getsockopt(conn, SOL_SOCKET, SO_SNDTIMEO, &old_timeout, &timeval_size) < 0) {
		goto err;
	}
	if (setsockopt(conn, SOL_SOCKET, SO_SNDTIMEO, timeout, sizeof(*timeout)) < 0) {
		goto err;
	}
	if (ws_ctube_socket_send_all(conn, response, strlen(response)) != 0) {
		goto err;
	}
	if (setsockopt(conn, SOL_SOCKET, SO_SNDTIMEO, &old_timeout, sizeof(old_timeout)) < 0) {
		goto err;
	}

	return 0;

err:
	if (WS_DEBUG) {
		perror("ws_ctube_ws_handshake()");
	}
	return -1;
}







#define WS_CTUBE_DEBUG 0
#define WS_CTUBE_BUFLEN 4096

typedef void (*cleanup_func)(void *);

static void _ws_ctube_cleanup_unlock_mutex(void *mutex)
{
	pthread_mutex_unlock((pthread_mutex_t *)mutex);
}

/** push a work item (start/stop connection) onto the FIFO connq */
static int ws_ctube_connq_push(struct ws_ctube *ctube, struct ws_ctube_conn_struct *conn, enum ws_ctube_qaction act)
{
	int retval = 0;
	struct ws_ctube_conn_qentry *qentry = (typeof(qentry))malloc(sizeof(*qentry));

	if (qentry == NULL) {
		retval = -1;
		goto out_noalloc;
	}
	pthread_cleanup_push(free, qentry);

	if (ws_ctube_conn_qentry_init(qentry, conn, act) != 0) {
		retval = -1;
		goto out_noinit;
	}
	pthread_cleanup_push((cleanup_func)ws_ctube_conn_qentry_destroy, qentry);

	ws_ctube_list_push_back(&ctube->connq, &qentry->lnode);
	pthread_mutex_lock(&ctube->connq_mutex);
	ctube->connq_pred = 1;
	pthread_mutex_unlock(&ctube->connq_mutex);
	pthread_cond_signal(&ctube->connq_cond);

	pthread_cleanup_pop(retval); /* conn_qentry_destory */
out_noinit:
	pthread_cleanup_pop(retval); /* free */
out_noalloc:
	return retval;
}

/** handles incoming data from client */
static void *ws_ctube_reader_main(void *arg)
{
	struct ws_ctube_conn_struct *conn = (struct ws_ctube_conn_struct *)arg;
	struct ws_ctube *ctube = conn->ctube;
	char buf[WS_CTUBE_BUFLEN];

	for (;;) {
		/* TODO: handle ping/pong */
		if (recv(conn->fd, buf, WS_CTUBE_BUFLEN, MSG_NOSIGNAL) < 1) {
			ws_ctube_connq_push(ctube, conn, WS_CTUBE_CONN_STOP);
			if (WS_CTUBE_DEBUG) {
				printf("ws_ctube_reader_main(): disconnected client\n");
				fflush(stdout);
			}
			return NULL;
		}
	}

	return NULL;
}

/** ensures data is released even if writer is cancelled */
static void _ws_ctube_cleanup_release_ws_ctube_data(void *arg)
{
	struct ws_ctube_data *ws_ctube_data = (struct ws_ctube_data *)arg;
	ws_ctube_ref_count_release(ws_ctube_data, refc, ws_ctube_data_free);
}

/** sends broadcast data to client */
static void *ws_ctube_writer_main(void *arg)
{
	struct ws_ctube_conn_struct *conn = (struct ws_ctube_conn_struct *)arg;
	struct ws_ctube *ctube = conn->ctube;
	struct ws_ctube_data *out_data = NULL;
	unsigned long out_data_id = 0;
	int send_retval;

	for (;;) {
		/* wait until new data is needed to be broadcast by checking data id */
		pthread_mutex_lock(&ctube->out_data_mutex);
		pthread_cleanup_push(_ws_ctube_cleanup_unlock_mutex, &ctube->out_data_mutex);
		while (out_data_id == ctube->out_data_id) {
			pthread_cond_wait(&ctube->out_data_cond, &ctube->out_data_mutex);
		}

		ws_ctube_ref_count_acquire(ctube->out_data, refc);
		out_data = ctube->out_data;
		out_data_id = ctube->out_data_id;

		pthread_cleanup_pop(0); /* _ws_ctube_cleanup_unlock_mutex */
		pthread_mutex_unlock(&ctube->out_data_mutex);

		/* broadcast data in a cancellable way */
		pthread_cleanup_push(_ws_ctube_cleanup_release_ws_ctube_data, out_data);
		send_retval = ws_ctube_ws_send(conn->fd, (char *)out_data->data, out_data->data_size);
		pthread_cleanup_pop(0); /* _ws_ctube_cleanup_release_ws_ctube_data */
		ws_ctube_ref_count_release(out_data, refc, ws_ctube_data_free);

		/* TODO: error handling of failed broadcast */
		if (send_retval != 0) {
			continue;
		}
	}

	return NULL;
}

static void _ws_ctube_cancel_reader(void *arg)
{
	int oldstate, statevar;
	pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &oldstate);

	struct ws_ctube_conn_struct *conn = (struct ws_ctube_conn_struct *)arg;
	pthread_cancel(conn->reader_tid);
	pthread_join(conn->reader_tid, NULL);

	pthread_setcancelstate(oldstate, &statevar);
}

static void _ws_ctube_cancel_writer(void *arg)
{
	int oldstate, statevar;
	pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &oldstate);

	struct ws_ctube_conn_struct *conn = (struct ws_ctube_conn_struct *)arg;
	pthread_cancel(conn->writer_tid);
	pthread_join(conn->writer_tid, NULL);

	pthread_setcancelstate(oldstate, &statevar);
}

/** start reader/writer threads for a client */
static int ws_ctube_conn_struct_start(struct ws_ctube_conn_struct *conn)
{
	int retval = 0;

	if (pthread_create(&conn->reader_tid, NULL, ws_ctube_reader_main, (void *)conn) != 0) {
		fprintf(stderr, "ws_ctube_conn_struct_start(): create reader failed\n");
		retval = -1;
		goto out_noreader;
	}
	pthread_cleanup_push(_ws_ctube_cancel_reader, conn);

	if (pthread_create(&conn->writer_tid, NULL, ws_ctube_writer_main, (void *)conn) != 0) {
		fprintf(stderr, "ws_ctube_conn_struct_start(): create writer failed\n");
		retval = -1;
		goto out_nowriter;
	}
	pthread_cleanup_push(_ws_ctube_cancel_writer, conn);

	pthread_cleanup_pop(retval); /* _ws_ctube_cancel_writer */
out_nowriter:
	pthread_cleanup_pop(retval); /* _ws_ctube_cancel_reader */
out_noreader:
	return retval;
}

/** cancels reader/writer threads for a client */
static void ws_ctube_conn_struct_stop(struct ws_ctube_conn_struct *conn)
{
	int oldstate, statevar;
	pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &oldstate);

	pthread_cancel(conn->reader_tid);
	pthread_cancel(conn->writer_tid);

	pthread_join(conn->reader_tid, NULL);
	pthread_join(conn->writer_tid, NULL);

	pthread_setcancelstate(oldstate, &statevar);
}

static void _ws_ctube_conn_list_add(struct ws_ctube_list *conn_list, struct ws_ctube_conn_struct *conn)
{
	ws_ctube_ref_count_acquire(conn, refc);
	ws_ctube_list_push_back(conn_list, &conn->lnode);
}

static void _ws_ctube_conn_list_remove(struct ws_ctube_list *conn_list, struct ws_ctube_conn_struct *conn)
{
	ws_ctube_list_unlink(conn_list, &conn->lnode);
	ws_ctube_ref_count_release(conn, refc, ws_ctube_conn_struct_free);
}

/** process work item from FIFO connq (start/stop connection) */
static void ws_ctube_handler_process_queue(struct ws_ctube_list *connq, struct ws_ctube_list *conn_list, int max_nclient)
{
	struct ws_ctube_conn_qentry *qentry;
	struct ws_ctube_list_node *node;
	struct ws_ctube_conn_struct *conn;

	while ((node = ws_ctube_list_pop_front(connq)) != NULL) {
		qentry = ws_ctube_container_of(node, typeof(*qentry), lnode);
		conn = qentry->conn;

		pthread_cleanup_push((cleanup_func)ws_ctube_conn_qentry_free, qentry);

		switch (qentry->act) {
		case WS_CTUBE_CONN_START:
			pthread_mutex_lock(&conn_list->mutex);

			/* refuse new connections if limit exceeded*/
			if (conn_list->len >= max_nclient) {
				pthread_mutex_unlock(&conn_list->mutex);
				fprintf(stderr, "ws_ctube_handler_process_queue(): max_nclient reached\n");
				fflush(stderr);
				break;
			} else {
				pthread_mutex_unlock(&conn_list->mutex);
			}

			/* do websocket handshake */
			if (ws_ctube_ws_handshake(conn->fd, &conn->ctube->timeout_val) == 0) {
				ws_ctube_conn_struct_start(conn);
				_ws_ctube_conn_list_add(conn_list, conn);
			}
			break;

		case WS_CTUBE_CONN_STOP:
			pthread_mutex_lock(&conn->stopping_mutex);
			/* prevent double stop */
			if (!conn->stopping) {
				conn->stopping = 1;
				pthread_mutex_unlock(&conn->stopping_mutex);

				_ws_ctube_conn_list_remove(conn_list, conn);
				ws_ctube_conn_struct_stop(conn);
			} else {
				pthread_mutex_unlock(&conn->stopping_mutex);
			}
			break;
		}

		pthread_cleanup_pop(1); /* ws_ctube_conn_qentry_free */
	}
}

/** stop all client connections and cleanup */
static void _ws_ctube_cleanup_conn_list(void *arg)
{
	int oldstate, statevar;
	pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &oldstate);

	struct ws_ctube_list *conn_list = (struct ws_ctube_list *)arg;
	struct ws_ctube_list_node *node;
	struct ws_ctube_conn_struct *conn;

	while ((node = ws_ctube_list_pop_front(conn_list)) != NULL) {
		conn = ws_ctube_container_of(node, typeof(*conn), lnode);

		pthread_mutex_lock(&conn->stopping_mutex);
		/* prevent double stop */
		if (!conn->stopping) {
			conn->stopping = 1;
			pthread_mutex_unlock(&conn->stopping_mutex);
			ws_ctube_conn_struct_stop(conn);
		} else {
			pthread_mutex_unlock(&conn->stopping_mutex);
		}

		ws_ctube_ref_count_release(conn, refc, ws_ctube_conn_struct_free);
	}

	pthread_setcancelstate(oldstate, &statevar);
}

/** connection handler thread: handles new clients via handshake or cleans up disconnected clients */
static void *ws_ctube_handler_main(void *arg)
{
	struct ws_ctube *ctube = (struct ws_ctube *)arg;

	struct ws_ctube_list conn_list;
	ws_ctube_list_init(&conn_list);
	pthread_cleanup_push(_ws_ctube_cleanup_conn_list, &conn_list);

	for (;;) {
		/* wait for work items in FIFO connq */
		pthread_mutex_lock(&ctube->connq_mutex);
		pthread_cleanup_push(_ws_ctube_cleanup_unlock_mutex, &ctube->connq_mutex);
		while (!ctube->connq_pred) {
			pthread_cond_wait(&ctube->connq_cond, &ctube->connq_mutex);
		}
		ctube->connq_pred = 0;
		pthread_mutex_unlock(&ctube->connq_mutex);
		pthread_cleanup_pop(0); /* _ws_ctube_cleanup_unlock_mutex */

		ws_ctube_handler_process_queue(&ctube->connq, &conn_list, ctube->max_nclient);
	}

	pthread_cleanup_pop(1);
	return NULL;
}

/* closes client socket on error */
static void _ws_ctube_cleanup_close_client_conn(void *arg)
{
	int oldstate, statevar;
	pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &oldstate);

	int *fd = (int *)arg;
	if (*fd >= 0) {
		close(*fd);
	}

	pthread_setcancelstate(oldstate, &statevar);
}

static int ws_ctube_server_accept_new_conn(struct ws_ctube *ctube, const int server_sock)
{
	int retval = 0;

	int conn_fd = accept(server_sock, NULL, NULL);
	pthread_cleanup_push(_ws_ctube_cleanup_close_client_conn, &conn_fd);

	/* create new conn_struct for client */
	struct ws_ctube_conn_struct *conn = (typeof(conn))malloc(sizeof(*conn));
	if (conn == NULL) {
		retval = -1;
		goto out_noalloc;
	}
	pthread_cleanup_push(free, conn);

	if (ws_ctube_conn_struct_init(conn, conn_fd, ctube) != 0) {
		retval = -1;
		goto out_noinit;
	}
	/* ws_ctube_conn_struct_destroy() closes conn_fd now; this prevents _ws_ctube_cleanup_close_client_conn() from closing it */
	conn_fd = -1;
	pthread_cleanup_push((cleanup_func)ws_ctube_conn_struct_destroy, conn);

	/* push work item to connection handler to do handshake and start reader/writer */
	if (ws_ctube_connq_push(ctube, conn, WS_CTUBE_CONN_START) != 0) {
		retval = -1;
		goto out_nopush;
	}

out_nopush:
	pthread_cleanup_pop(retval); /* ws_ctube_conn_struct_destroy */
out_noinit:
	pthread_cleanup_pop(retval); /* free */
out_noalloc:
	pthread_cleanup_pop(retval); /* _ws_ctube_cleanup_close_client_conn */
	return retval;
}

static void ws_ctube_serve_forever(struct ws_ctube *ctube)
{
	const int server_sock = ctube->server_sock;

	for (;;) {
		if (ws_ctube_server_accept_new_conn(ctube, server_sock) != 0) {
			fprintf(stderr, "ws_ctube_serve_forever(): error\n");
			fflush(stderr);
		}
	}
}

/** alert main thread if the server fails to init by setting flag */
static void _ws_ctube_server_init_fail(void *arg)
{
	int oldstate, statevar;
	pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &oldstate);

	struct ws_ctube *ctube = (struct ws_ctube *)arg;
	pthread_mutex_lock(&ctube->server_init_mutex);
	ctube->server_inited = -1;
	pthread_mutex_unlock(&ctube->server_init_mutex);
	pthread_cond_signal(&ctube->server_init_cond);

	pthread_setcancelstate(oldstate, &statevar);
}

static void _ws_ctube_close_server_sock(void *arg)
{
	int oldstate, statevar;
	pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &oldstate);

	struct ws_ctube *ctube = (struct ws_ctube *)arg;
	close(ctube->server_sock);

	pthread_setcancelstate(oldstate, &statevar);
}

static void *ws_ctube_server_main(void *arg)
{
	struct ws_ctube *ctube = (struct ws_ctube *)arg;
	pthread_cleanup_push(_ws_ctube_server_init_fail, ctube);

	/* create server socket */
	int server_sock = socket(AF_INET, SOCK_STREAM, 0);
	if (server_sock < 0) {
		perror("ws_ctube_server_main()");
		goto out_nosock;
	}
	ctube->server_sock = server_sock;
	pthread_cleanup_push(_ws_ctube_close_server_sock, ctube);

	/* allow reuse */
#ifdef __linux__
	int optname = SO_REUSEADDR | SO_REUSEPORT;
#else
	int optname = SO_REUSEADDR;
#endif
	int yes = 1;
	if (setsockopt(server_sock, SOL_SOCKET, optname, &yes, sizeof(yes)) < 0) {
		perror("ws_ctube_server_main()");
		goto out_err;
	}

	/* set server socket address/port */
	if (ws_ctube_bind_server(server_sock, ctube->port) < 0) {
		perror("ws_ctube_server_main()");
		goto out_err;
	}

	/* set listening */
	if (listen(server_sock, ctube->max_nclient) < 0) {
		perror("ws_ctube_server_main()");
		goto out_err;
	}

	/* success: alert main thread by setting flag */
	pthread_mutex_lock(&ctube->server_init_mutex);
	ctube->server_inited = 1;
	pthread_mutex_unlock(&ctube->server_init_mutex);
	pthread_cond_signal(&ctube->server_init_cond);
	ws_ctube_serve_forever(ctube);

	/* code doesn't get here unless error */
out_err:
	pthread_cleanup_pop(1); /* _ws_ctube_close_server_sock */
out_nosock:
	pthread_cleanup_pop(1); /* _ws_ctube_server_init_fail */
	return NULL;
}

static void _ws_ctube_cancel_handler(void *arg)
{
	int oldstate, statevar;
	pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &oldstate);

	struct ws_ctube *ctube = (struct ws_ctube *)arg;
	pthread_cancel(ctube->handler_tid);
	pthread_join(ctube->handler_tid, NULL);

	pthread_setcancelstate(oldstate, &statevar);
}

static void _ws_ctube_cancel_server(void *arg)
{
	int oldstate, statevar;
	pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &oldstate);

	struct ws_ctube *ctube = (struct ws_ctube *)arg;
	pthread_cancel(ctube->server_tid);
	pthread_join(ctube->server_tid, NULL);

	pthread_setcancelstate(oldstate, &statevar);
}

/* start connection handler and server threads */
static int ws_ctube_start(struct ws_ctube *ctube)
{
	int retval = 0;

	if (pthread_create(&ctube->handler_tid, NULL, ws_ctube_handler_main, (void *)ctube) != 0) {
		fprintf(stderr, "ws_ctube_start(): create handler failed\n");
		retval = -1;
		goto out_nohandler;
	}
	pthread_cleanup_push(_ws_ctube_cancel_handler, ctube);

	if (pthread_create(&ctube->server_tid, NULL, ws_ctube_server_main, (void *)ctube) != 0) {
		fprintf(stderr, "ws_ctube_start(): create server failed\n");
		retval = -1;
		goto out_noserver;
	}
	pthread_cleanup_push(_ws_ctube_cancel_server, ctube);

	/* wait for server thread to report success/failure to start */
	pthread_mutex_lock(&ctube->server_init_mutex);
	pthread_cleanup_push(_ws_ctube_cleanup_unlock_mutex, &ctube->server_init_mutex);
	if (ctube->timeout_spec.tv_nsec > 0 || ctube->timeout_spec.tv_sec > 0) {
		while (!ctube->server_inited) {
			pthread_cond_timedwait(&ctube->server_init_cond, &ctube->server_init_mutex, &ctube->timeout_spec);
		}
	} else {
		while (!ctube->server_inited) {
			pthread_cond_wait(&ctube->server_init_cond, &ctube->server_init_mutex);
		}
	}
	if (ctube->server_inited <= 0) {
		fprintf(stderr, "ws_ctube_start(): server failed to init\n");
		retval = -1;
		goto out_noinit;
	}
	pthread_mutex_unlock(&ctube->server_init_mutex);

out_noinit:
	pthread_cleanup_pop(retval); /* _ws_ctube_cleanup_unlock_mutex */
	pthread_cleanup_pop(retval); /* _ws_ctube_cancel_server */
out_noserver:
	pthread_cleanup_pop(retval); /* _ws_ctube_cancel_handler */
out_nohandler:
	return retval;
}

/** stop connection handler and server threads */
static void ws_ctube_stop(struct ws_ctube *ctube)
{
	int oldstate, statevar;
	pthread_setcancelstate(PTHREAD_CANCEL_DISABLE, &oldstate);

	pthread_cancel(ctube->handler_tid);
	pthread_cancel(ctube->server_tid);

	pthread_join(ctube->handler_tid, NULL);
	pthread_join(ctube->server_tid, NULL);

	pthread_setcancelstate(oldstate, &statevar);
}

struct ws_ctube* ws_ctube_open(
	int port,
	int max_nclient,
	int timeout_ms,
	double max_broadcast_fps)
{
	int err = 0;
	struct ws_ctube *ctube;

	/* input sanity checks */
	if (port < 1) {
		fprintf(stderr, "ws_ctube_open(): invalid port\n");
		fflush(stderr);
		err = -1;
		goto out_noalloc;
	}
	if (max_nclient < 1) {
		fprintf(stderr, "ws_ctube_open(): invalid max_nclient\n");
		fflush(stderr);
		err = -1;
		goto out_noalloc;
	}
	if (timeout_ms < 0) {
		fprintf(stderr, "ws_ctube_open(): invalid timeout_ms\n");
		fflush(stderr);
		err = -1;
		goto out_noalloc;
	}
	if (max_broadcast_fps < 0) {
		fprintf(stderr, "ws_ctube_open(): invalid max_broadcast_fps\n");
		fflush(stderr);
		err = -1;
		goto out_noalloc;
	}

	ctube = (typeof(ctube))malloc(sizeof(*ctube));
	if (ctube == NULL) {
		err = -1;
		goto out_noalloc;
	}
	pthread_cleanup_push((cleanup_func)free, ctube);

	if (ws_ctube_init(ctube, port, max_nclient, timeout_ms, max_broadcast_fps) != 0) {
		err = -1;
		goto out_noinit;
	}
	pthread_cleanup_push((cleanup_func)ws_ctube_destroy, ctube);

	if (ws_ctube_start(ctube) != 0) {
		err = -1;
		goto out_nostart;
	}

out_nostart:
	pthread_cleanup_pop(err);
out_noinit:
	pthread_cleanup_pop(err);
out_noalloc:
	if (err) {
		return NULL;
	} else {
		return ctube;
	}
}

void ws_ctube_close(struct ws_ctube *ctube)
{
	if (ws_ctube_unlikely(ctube == NULL)) {
		fprintf(stderr, "ws_ctube_close(): error: ctube is NULL\n");
		fflush(stderr);
		return;
	}

	ws_ctube_stop(ctube);
	ws_ctube_destroy(ctube);
	free(ctube);
}

int ws_ctube_broadcast(struct ws_ctube *ctube, const void *data, size_t data_size)
{
	if (ws_ctube_unlikely(ctube == NULL)) {
		fprintf(stderr, "ws_ctube_broadcast(): error: ctube is NULL\n");
		fflush(stderr);
		return -1;
	}
	if (ws_ctube_unlikely(data == NULL)) {
		fprintf(stderr, "ws_ctube_broadcast(): error: data is NULL\n");
		fflush(stderr);
		return -1;
	}
	if (ws_ctube_unlikely(data_size == 0)) {
		fprintf(stderr, "ws_ctube_broadcast(): error: data_size is 0\n");
		fflush(stderr);
		return -1;
	}

	int retval = 0;
	if (pthread_mutex_trylock(&ctube->out_data_mutex) != 0) {
		retval = -1;
		goto out_nolock;
	}
	pthread_cleanup_push(_ws_ctube_cleanup_unlock_mutex, &ctube->out_data_mutex);

	/* rate limit broadcasting if set */
	struct timespec cur_time;
	const double max_bcast_fps = ctube->max_bcast_fps;
	if (max_bcast_fps > 0) {

#ifdef CLOCK_MONOTONIC
		if (ws_ctube_unlikely(clock_gettime(CLOCK_MONOTONIC, &cur_time) != 0)) {
			clock_gettime(CLOCK_REALTIME, &cur_time);
		}
#else
		clock_gettime(CLOCK_REALTIME, &cur_time);
#endif /* CLOCK_MONOTONIC */

		double dt = (cur_time.tv_sec - ctube->prev_bcast_time.tv_sec) +
			1e-9 * (cur_time.tv_nsec - ctube->prev_bcast_time.tv_nsec);

		if (dt < 1.0 / ctube->max_bcast_fps) {
			retval = -1;
			goto out_ratelim;
		}
	}

	/* release old out_data if held */
	if (ctube->out_data != NULL) {
		ws_ctube_ref_count_release(ctube->out_data, refc, ws_ctube_data_free);
	}

	/* alloc new out_data */
	ctube->out_data = (typeof(ctube->out_data))malloc(sizeof(*ctube->out_data));
	if (ws_ctube_unlikely(ctube->out_data == NULL)) {
		retval = -1;
		goto out_nodata;
	}
	pthread_cleanup_push(free, ctube->out_data);

	/* init and memcpy into out_data */
	if (ws_ctube_unlikely(ws_ctube_data_init(ctube->out_data, data, data_size) != 0)) {
		retval = -1;
		goto out_noinit;
	}
	ws_ctube_ref_count_acquire(ctube->out_data, refc);
	ctube->out_data_id++; /* unique id for out_data */

	/* record broadcast time for rate-limiting next time */
	if (max_bcast_fps > 0) {
		ctube->prev_bcast_time = cur_time;
	}

	pthread_mutex_unlock(&ctube->out_data_mutex);
	pthread_cond_broadcast(&ctube->out_data_cond);

out_noinit:
	pthread_cleanup_pop(retval); /* free */
out_nodata:
out_ratelim:
	pthread_cleanup_pop(retval); /* _ws_ctube_cleanup_unlock_mutex */
out_nolock:
	return retval;
}


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */


#endif /* WS_CTUBE_H */

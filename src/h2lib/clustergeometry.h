#ifndef CLUSTERGEOMETRY_H
#define CLUSTERGEOMETRY_H

typedef struct _clustergeometry clustergeometry;
typedef clustergeometry* pclustergeometry;

#include "cluster.h"
#include "h2settings.h"
struct _clustergeometry {
    /** @brief Spatial dimension. */
    int dim;

    /** @brief Number of indices.*/
    int nidx;

    /**  @brief Characteristic points for the indices.*/
    double **x;

    /** @brief Minimal coordinates for the support bounding boxes for the indices.*/
    double **smin;

    /** @brief Maximal coordinates for the support bounding boxes for the indices.*/
    double **smax;

    /** @brief Weights for the indices. */
    double *w;

    /** @brief internal fields used to build the @ref cluster tree from the
    clustergeometry object.*/
    double *hmin;

    /** @brief internal fields used to build the @ref cluster tree from the
    clustergeometry object.*/
    double *hmax;

    /** @brief internal fields used to build the @ref cluster tree from the
    clustergeometry object.*/
    double *buf;
};

/* ------------------------------------------------------------
 Constructors and destructors
 ------------------------------------------------------------ */

/**
 * @brief Create a new @ref clustergeometry object.
 *
 * Allocates storage for the object and the internal fields.
 *
 * @remark Should always be matched by a call to @ref del_clustergeometry.
 *
 * @param dim Spatial dimension of the domain.
 * @param nidx Number of characteristic points in the domain.
 * @return Returns the newly created @ref clustergeometry object.
 */
pclustergeometry
new_clustergeometry(int dim, int nidx);

/**
 * @brief Delete a @ref clustergeometry object.
 *
 * Releases the storage corresponding to this object.
 *
 *@param cf Object to be deleted. */
void
del_clustergeometry(pclustergeometry cf);

/* ------------------------------------------------------------
 Auxiliary routines
 ------------------------------------------------------------ */

/**
 * @brief Update an adaptive bounding box for an index set.
 *
 *  Computes an adaptive bounding box for an index set and stores the result
 *  in the internal fields <tt>hmin</tt> and <tt>hmax</tt> of the
 *  @ref clustergeometry structure.
 *  Already existing entries are overwritten.
 *
 *  @param cf clustergeometry object, where the bounding box is stored.
 *  @param size Number of indices.
 *  @param idx Index set. */
void
update_bbox_for_clustergeometry(pclustergeometry cf, int size, int *idx);

/**
 * @brief Update a bounding box for the support of a cluster.
 *
 *  Computes a bounding box for the support of a cluster <tt> t </tt>using only
 *  the fields <tt>smin</tt> and <tt>smax</tt> and the index set <tt>idx</tt>
 *  of the @ref clustergeometry structure.
 *  Already existing entries are overwritten.
 *
 *  @param cf Clustergeometry object with geometrical information.
 *  @param t In this @ref cluster tree the bounding boxes are updated. */
void
update_sbox_for_cluster(pclustergeometry cf, pcluster t);

/**
 * @brief Updates the bounding boxes in a @ref cluster tree object using
 *  only its sons.
 *
 *  Computes the coordinates of the support bounding boxes of a cluster using
 *  only the bounding boxes of its sons.
 *  Already existing values are overwritten.
 *
 *  @param t In this cluster the bounding boxes are updated. */
void
update_bbox_for_cluster(pcluster t);

#endif

/** @}*/

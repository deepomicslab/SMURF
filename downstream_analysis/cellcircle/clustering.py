# -*- coding: utf-8 -*-
"""
Created on Sat May 14 17:05:57 2022

@author: WANG Bingchen
"""

from sklearn.cluster import SpectralClustering
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn import metrics


def spectral_clustering(X, k):
    label = SpectralClustering(
        n_clusters=k,
        affinity="nearest_neighbors",
        n_neighbors=40,
        random_state=0
    ).fit_predict(X)
    return label


def kmeans_clustering(X, k):
    label = KMeans(n_clusters=k).fit_predict(X)
    return label


def hierarchical_clustering(X, k):
    label = AgglomerativeClustering(
        n_clusters=k,
        affinity='euclidean',
        linkage='ward',
        ).fit_predict(X)
    return label


def get_metric(label, pred_label):
    return {
        'ARI': metrics.adjusted_rand_score(label, pred_label),
        'AMI': metrics.adjusted_mutual_info_score(label, pred_label),
        'NMI': metrics.normalized_mutual_info_score(label, pred_label),
        'Homogeneity': metrics.homogeneity_score(label, pred_label),
        'Completeness': metrics.completeness_score(label, pred_label),
        'V-measure': metrics.v_measure_score(label, pred_label),

    }


def get_metric_without_label(X, label):
    return {
        'SW': metrics.silhouette_score(X, label, metric='cosine'),
        # 'CHI': metrics.calinski_harabasz_score(X, label),
        # 'DBI': metrics.davies_bouldin_score(X, label),
    }
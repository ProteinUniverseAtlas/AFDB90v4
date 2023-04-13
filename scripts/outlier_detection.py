import itertools
import pickle
from multiprocessing import Pool
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from gensim.models import FastText
from sklearn.ensemble import IsolationForest
from tqdm import tqdm
from collections import defaultdict, Counter
import pickle

global VECTORIZER

def documents_to_keys_corpus(documents):
    keys_corpus = (line.strip().split("\t") for corpus_file in documents for line in tqdm(open(corpus_file)) if
                   len(line.strip().split("\t")) == 2)
    keys, corpus = itertools.tee(keys_corpus)
    keys = [k[0] for k in keys]
    corpus = (k[1] for k in corpus)
    return keys, corpus


def vectorizer_from_documents(documents, fasttext_model_file, corpus_file, num_jobs=100):
    keys, corpus = documents_to_keys_corpus(documents)
    with open(corpus_file, "w") as f:
        for shapemers in corpus:
            f.write(" ".join(f'{int(s):010b}' for s in shapemers.split()) + "\n")
    vectorizer = FastText(vector_size=1024, window=16, min_count=1, workers=num_jobs)
    vectorizer.build_vocab(corpus_file=corpus_file)
    vectorizer.train(
        corpus_file=corpus_file, epochs=vectorizer.epochs,
        total_examples=vectorizer.corpus_count, total_words=vectorizer.corpus_total_words,
    )
    vectorizer.save(str(fasttext_model_file))
    
    
def get_sentence_vector(document):
    return VECTORIZER.get_sentence_vector([f'{int(s):010b}' for s in document.split()])

def vectorize(corpus, num_jobs, total):
    with Pool(processes=num_jobs) as pool:
        matrix = list(tqdm(pool.imap(get_sentence_vector,
                                     corpus),
                           total=total))
    return np.array(matrix)

def make_outlier_figure():
    data_file = pnd.read_csv("data/data.csv", index_col=0)
    to_darkness = dict(zip(data_file["AF2_longest_best70"], data_file["FULL_noDUF"]))
    to_length = dict(zip(data_file["AF2_longest_best70"], data_file["AF2_longest_best70_len"]))
    to_diversity = dict(zip(data_file["AF2_longest_best70"], data_file["diversity_fraction"]))
    to_outlier = {}
    with open("data/afdb90_outliers.txt") as f:
        for line in tqdm(f):
            key, score = line.strip().split("\t")
            to_outlier[key] = float(score)
    outlier_color = "#634248"
    inlier_color = "#758a5e"
    background_color = '#F2F2F2'
    purple = '#57257F'
    inlier_threshold = 0.115
    fig, ax = plt.subplots(2, 3, figsize=(10, 10), sharey=True)
    for i in range(2):
        for j in range(3):
            ax[i][j].set_facecolor(background_color)

    nbins = np.arange(0, 105, 5)
    n = list(range(0,21,5))
    l = list(range(0,101,25))
    heights, _ = np.histogram([to_darkness[k[:-3]] for k in to_outlier if to_outlier[k] < 0], bins = nbins)
    heights = heights * 100/sum(heights)
    x = list(range(len(heights)))
    y = list(heights)
    ax[0][0].bar(x, y, 1, align="edge", 
             label="Outliers", color=outlier_color, edgecolor="black");
    ax[0][0].set_xlabel("Functional brightness (%)")
    ax[0][0].set_xticks(n, l)

    heights, _ = np.histogram([to_darkness[k[:-3]] for k in to_outlier if to_outlier[k] > inlier_threshold], bins = nbins)
    heights = heights * 100/sum(heights)
    x = list(range(len(heights)))
    y = list(heights)
    ax[1][0].bar(x, y, 1, align="edge", 
             label="Inliers", color=inlier_color, edgecolor="black");
    ax[1][0].set_xlabel("Functional brightness (%)")
    ax[1][0].set_xticks(n, l)


    nbins = np.arange(0, 1.1, 0.1)
    n = [(i, f"{nbins[i]:.1f}") for i in range(len(nbins)) if i%2==0]
    heights, _ = np.histogram([to_diversity[k[:-3]] for k in to_outlier if to_outlier[k] < 0], bins = nbins)
    heights = heights * 100/sum(heights)
    x = list(range(len(heights)))
    y = list(heights)
    ax[0][1].bar(x, y, 1, align="edge", 
             label="Outliers", color=outlier_color, edgecolor="black");
    ax[0][1].set_xlabel("$Diversity=\\frac{|\{shape-mers\}|}{|shape-mers|}$")
    ax[0][1].set_xticks([x[0] for x in n], [x[1] for x in n])

    heights, _ = np.histogram([to_diversity[k[:-3]] for k in to_outlier if to_outlier[k] > inlier_threshold], bins = nbins)
    heights = heights * 100/sum(heights)
    x = list(range(len(heights)))
    y = list(heights)
    ax[1][1].bar(x, y, 1, align="edge", 
             label="Inliers", color=inlier_color, edgecolor="black");
    ax[1][1].set_xlabel("$Diversity=\\frac{|\{shape-mers\}|}{|shape-mers|}$")
    ax[1][1].set_xticks([x[0] for x in n], [x[1] for x in n])

    ax[0][0].set_title("A", fontsize=15, loc="left")
    ax[0][1].set_title("Outliers", fontsize=15)
    ax[1][0].set_title("B", fontsize=15, loc="left")
    ax[1][1].set_title("Inliers", fontsize=15)

    nbins = np.arange(0, 2250, 250)
    n = [(i, nbins[i]) for i in range(len(nbins)) if i%2==0]
    heights, _ = np.histogram([to_length[k[:-3]] for k in to_outlier if to_outlier[k] < 0], bins = nbins)
    heights = heights * 100/sum(heights)
    x = list(range(len(heights)))
    y = list(heights)
    ax[0][2].bar(x, y, 1, align="edge", 
             label="Outliers", color=outlier_color, edgecolor="black");
    ax[0][2].set_xlabel("Protein length")
    ax[0][2].set_xticks([x[0] for x in n], [x[1] for x in n])

    heights, _ = np.histogram([to_length[k[:-3]] for k in to_outlier if to_outlier[k] > inlier_threshold], bins = nbins)
    heights = heights * 100/sum(heights)
    x = list(range(len(heights)))
    y = list(heights)
    ax[1][2].bar(x, y, 1, align="edge", 
             label="Inliers", color=inlier_color, edgecolor="black");
    ax[1][2].set_xlabel("Protein length")
    ax[1][2].set_xticks([x[0] for x in n], [x[1] for x in n])

    ax[0][0].set_ylabel("% of outliers")
    ax[1][0].set_ylabel("% of inliers")
    plt.subplots_adjust(wspace=0.1, hspace=0.3)
    plt.yticks(np.arange(0, 110, 10))
    plt.savefig("outliers.pdf", dpi=2000)
    

def main():
    num_jobs = 24
    vectorizer_from_documents(["data/shapemers_proteinnet/proteinnet_shapemers.txt"], 
                          "data/fasttext.mdl", 
                          "data/fasttext_corpus.txt", num_jobs=num_jobs)
    VECTORIZER = FastText.load("data/fasttext.mdl", mmap='r').wv
    outlier_detector = IsolationForest(n_jobs=num_jobs, n_estimators=1000,
                                       max_features=0.5, contamination=0.1,
                                       verbose=10)
    keys, corpus = documents_to_keys_corpus(["data/shapemers_proteinnet/proteinnet_shapemers.txt"])
    matrix = vectorize(corpus, num_jobs, len(keys))
    outlier_detector.fit(matrix)
    with open("data/isolation_forest.pkl", "wb") as f:
        pickle.dump(outlier_detector, f)
    output_folder = Path("data/outlier_results")
    output_folder.mkdir(exist_ok=True)
    for filename in Path("data/afdb50_file_split/").glob("*.txt"):
        with open(filename, "r") as f:
            documents = [Path("data/shapemers_afdb50/") / line.strip() for line in f]
        keys, corpus = documents_to_keys_corpus(documents)
        print("Vectorizing")
        matrix = np.array(
            [vectorizer.get_sentence_vector([f'{int(s):010b}' for s in document.split()]) for document in tqdm(corpus)])
        print("Predicting")
        scores = outlier_detector.decision_function(matrix)
        print("Saving")
        with open(output_folder / f"{filename.stem}.txt", "w") as f:
            for key, score in zip(keys, scores):
                f.write(f"{key}\t{score:.3f}\n")
    
    
if __name__ == "__main__":
    main()
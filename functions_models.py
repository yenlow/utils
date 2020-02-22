import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)


def sk_mod_summary(fit, features):
    """
    Returns coef_ of sklearn's fit() in a pandas series for easy viewing
    Args:
        fit: model class from calling fit()
        features: list of features

    Returns:

    """
    return pd.Series(index=features, data=fit.coef_[0].tolist())



def sk_pipeline_step_names(pl, match='classifier'):
    """
    Returns step names in pipeline and finds the number of step matching a provided name (e.g. 'classifier' as default)
    Args:
        pl (object): sklearn.pipeline.Pipeline class
        match (str): step name to be located (e.g. 'classifier' as default)

    Returns:
        stepnames as list, position of match in stepnames list

    Example:
        >>>> sk_pipeline_step_names(pl, 'classifier')
        (['preprocessor', 'classifier'], 1)

    """
    step_names = [name for (name, obj) in pl.steps]
    if match:
        pos = step_names.index(match)
    else:
        pos = None

    return step_names, pos


def bestLambda_lrcv(lrcv, plot=True, title='Max score from CV'):
    """
    Returns best hyperparameter, score type and pandas dataframe of scores for hyperparameters (columns) per fold (rows)
    Also plots scores vs hyperparameters with each line representing sample from each CV fold
    Args:
        lrcv: object from LogisticRegressionCV()
    Returns:
        best_C (float): best hyperparameter
        score_type (str):  optimization metric (accessed by lrcv.scoring)
        scores_df (pd.DataFrame): scores for hyperparameters (columns) per fold (rows)

    """
    best_C = lrcv.C_[0]
    score_type = lrcv.scoring
    scores_df = pd.DataFrame(lrcv.scores_[1], columns=lrcv.Cs_)

    if plot:
        scores_df.T.plot.line(logx=True, title=title)
        plt.xlabel('C or 1/lambda')
        plt.ylabel(lrcv.scoring)
        plt.axvline(x=best_C, linewidth=2, linestyle='dotted', c='r')
        plt.text(best_C, 0.9, str(best_C)[0:6], fontsize=12)
        plt.legend().set_bbox_to_anchor((10e3, lrcv.scores_[1].min()))
        plt.show()
    return best_C, score_type, scores_df


# https://stackoverflow.com/questions/25122999/scikit-learn-how-to-check-coefficients-significance
# Wald's approximation
# TODO: illconditioning
def logit_pvalue(model, x):
    """ Calculate z-scores for scikit-learn LogisticRegression.
    parameters:
        model: fitted sklearn.linear_model.LogisticRegression with intercept and large C
        x:     matrix on which the model was fit
    This function uses asymtptics for maximum likelihood estimates.
    """
    p = model.predict_proba(x)
    n = len(p)
    m = len(model.coef_[0]) + 1
    coefs = np.concatenate([model.intercept_, model.coef_[0]])
    x_full = np.matrix(np.insert(np.array(x), 0, 1, axis = 1))
    ans = np.zeros((m, m))
    for i in range(n):
        ans = ans + np.dot(np.transpose(x_full[i, :]), x_full[i, :]) * p[i,1] * p[i, 0]
    vcov = np.linalg.inv(np.matrix(ans))
    se = np.sqrt(np.diag(vcov))
    t =  coefs/se
    p = (1 - norm.cdf(abs(t))) * 2
    return p
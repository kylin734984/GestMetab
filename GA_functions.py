import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import cross_val_score, KFold
from deap import base, creator, tools, algorithms
import multiprocessing
from sklearn.exceptions import ConvergenceWarning
from sklearn.utils._testing import ignore_warnings


@ignore_warnings(category=ConvergenceWarning)
def evalOneMax(individual, X_train, y_train, nvar):
    selected_features = [index for index in range(len(individual)) if individual[index] == 1]
    if len(selected_features) != nvar:
        return 1E200,  # 返回一个大错误，以防止选择的特征数量不等于指定数量
    X_selected = X_train[:, selected_features]
    model = MLPRegressor(hidden_layer_sizes=(32,), activation='relu', solver='adam', max_iter=100, random_state=123)
    kf = KFold(n_splits=3, shuffle=True, random_state=123)
    scores = cross_val_score(model, X_selected, y_train, cv=kf, scoring='neg_mean_squared_error')
    print('evaluate...')
    return -1 * np.mean(scores),


def mutate_and_fix(individual, indpb):
    tools.mutFlipBit(individual, indpb)
    while sum(individual) > 20:
        idx = np.random.choice(np.where(np.array(individual) == 1)[0])
        individual[idx] = 0
    while sum(individual) < 20:
        idx = np.random.choice(np.where(np.array(individual) == 0)[0])
        individual[idx] = 1
    return individual,


def ga_nvar_rf_reg(X_train, y_train, workers=1, initpop=50, nvar=20, cxpb=0.5, mutpb=0.2, ngen=50):
    # 遗传算法进行特征选择
    creator.create("FitnessMin", base.Fitness, weights=(-1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMin)

    toolbox = base.Toolbox()
    toolbox.register("attr_bool", np.random.choice, [0, 1])
    toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_bool, n=len(X_train[0]))
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    toolbox.register("evaluate", evalOneMax, X_train=X_train, y_train=y_train, nvar=nvar)
    toolbox.register("mate", tools.cxTwoPoint)
    toolbox.register("mutate", mutate_and_fix, indpb=0.05)
    toolbox.register("select", tools.selTournament, tournsize=3)

    # 使用多处理
    pool = multiprocessing.Pool(processes=workers)
    toolbox.register("map", pool.map)

    # 初始化种群
    population = toolbox.population(n=initpop)
    for individual in population:
        while sum(individual) > nvar:
            idx = np.random.choice(np.where(np.array(individual) == 1)[0])
            individual[idx] = 0
        while sum(individual) < nvar:
            idx = np.random.choice(np.where(np.array(individual) == 0)[0])
            individual[idx] = 1

    result, logbook = algorithms.eaSimple(population, toolbox, cxpb, mutpb, ngen, verbose=False)
    pool.close()
    pool.join()

    # 选择最好的特征
    best_individual = tools.selBest(population, k=1)[0]
    selected_features = [index for index in range(len(best_individual)) if best_individual[index] == 1]
    return selected_features

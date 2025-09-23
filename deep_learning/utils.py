import os
import yaml
import h5py
import pickle
import torch
import torch.nn as nn
import matplotlib.pyplot as plt
from datetime import datetime

from sklearn.metrics import (
    confusion_matrix,
    accuracy_score,
    precision_score,
    recall_score,
    f1_score)

from engine_setup import initialize_metrics


def load_model_input_data(data_directory_path: str, parameters: dict) -> tuple:
    data_file_path = os.path.join(
        data_directory_path,
        f"{parameters['data_version']}.hdf5")
    
    with h5py.File(data_file_path, "r") as file:
        images = file["X"][:]
        labels = file["Y"][:]

    images = images.transpose(0, 4, 1, 2, 3)

    return images, labels


def generate_directory_path(target_path: str, prefix: str = None) -> str:
    timestamp = get_timestamp()
    directory_name = f"{prefix}_{timestamp}"
    return os.path.join(target_path, directory_name)


def ensure_directory_exists(target_path: str) -> str:
    os.makedirs(target_path, exist_ok=True)


def get_timestamp():
    return datetime.now().strftime("%Y-%m-%d_%H-%M")


def save_config_params(target_path: str, parameters: dict, file_name: str) -> None:
    with open(os.path.join(target_path, file_name), "w") as file:
        yaml.dump(parameters, file, default_flow_style=False)


def log_params_saved(file_name: str) -> None:
    print(f"Configuration file '{file_name}' has been successfully saved.")


def log_saved_file_path(target_path: str) -> None:
    print(f"Path to saved file(s): {target_path}")


def get_device_from_model(model):
    return next(model.parameters()).device


def transfer_to_device(tensor, device):
    return tensor.to(device)


def transfer_data_to_device(features, targets, device):
    features = transfer_to_device(features, device)
    targets = transfer_to_device(targets, device)

    return features, targets


def log_training_start():
    print("Training has started...\n")


def log_epoch_results(epoch, train_metrics, valid_metrics):
    print(f"\nEpoch: {epoch}")
    print(f"Training loss: {train_metrics[0]:.4f} | Validation loss: {valid_metrics[0]:.4f}")
    print(f"Training accuracy: {train_metrics[1]:.4f} | Validation accuracy: {valid_metrics[1]:.4f}\n")


def plot_training_curves(results: dict, metric: str) -> plt.Figure:
    num_epochs = len(results[metric]["train"])
    epoch_range = range(1, num_epochs + 1)
    
    fig, ax = plt.subplots()
    ax.plot(epoch_range, results[metric]["train"], label=f"Training {metric.capitalize()} Curve")
    ax.plot(epoch_range, results[metric]["valid"], label=f"Validation {metric.capitalize()} Curve")
    ax.set_title(f"{metric.capitalize()} Curves for Training and Validation Data")
    ax.set_xlabel("Epochs")
    ax.set_ylabel(f"{metric.capitalize()}")
    ax.grid(True)
    ax.legend()

    return fig


def prepare_directory(target_path: str, directory_name: str) -> str:
    directory_path = generate_directory_path(target_path, directory_name)
    ensure_directory_exists(directory_path) 

    return directory_path


def handle_training_artifacts_saving(
        model: nn.Module,
        checkpoints: dict,
        training_parameters: dict,
        evaluation_report: dict,
        target_path: str,
) -> None:
    directory_path = prepare_directory(
        target_path, training_parameters["model_name"])
    

    save_model(directory_path, model)
    log_model_saved()

    save_checkpoints(directory_path, checkpoints)
    log_checkpoints_saved()

    training_metrics = checkpoints["training_metrics"]
    loss_figure = plot_training_curves(training_metrics, metric="loss")
    save_plot(directory_path, loss_figure, file_name="loss.png")

    accuracy_figure = plot_training_curves(training_metrics, metric="accuracy")
    save_plot(directory_path, accuracy_figure, file_name="accuracy.png")

    print_training_summary(training_metrics)
    save_training_summary(directory_path, training_metrics)

    save_evaluation_report(directory_path, evaluation_report)

    config_file_name = "model_parameters.txt"
    save_config_params(directory_path, training_parameters, config_file_name)
    log_params_saved(config_file_name)

    log_saved_file_path(directory_path)


def save_model(target_dir, model, file_name="model.pth"):
    torch.save(model.state_dict(), os.path.join(target_dir, file_name))


def log_model_saved():
    print("Model has been successfully saved.")

    
def create_checkpoints(optimizer, scheduler, metrics, best_step) -> dict:
    return {
        "optimizer": optimizer,
        "lr_scheduler": scheduler,
        "training_metrics": metrics,
        "best_training_step": best_step
    }


def save_checkpoints(target_dir, checkpoints, file_name="my_checkpoints.pth"):
    torch.save(checkpoints, os.path.join(target_dir, file_name))


def log_checkpoints_saved():
    print("Checkpoints have been successfully saved.")


def save_plot(target_dir, figure, file_name="plot.png"):
    figure.savefig(os.path.join(target_dir, file_name))


def log_plot_saved():
    print("Plot has been successfully saved.")


def initialize_results_tracker() -> dict:
    return {
        "loss": {
            "train": [],
            "valid": [],
        },
        "accuracy": {
            "train": [],
            "valid": [],
        },
    }


def update_results_tracker(
    results: dict, train_metrics: tuple, valid_metrics: tuple) -> None:
    results["loss"]["train"].append(train_metrics[0])
    results["loss"]["valid"].append(valid_metrics[0])
    results["accuracy"]["train"].append(train_metrics[1])
    results["accuracy"]["valid"].append(valid_metrics[1])


def is_stopper_triggered(stopper, valid_loss: float) -> bool:
    return stopper and stopper.stop(valid_loss)


def log_early_stopping(stopper_name: str) -> None:
    print(f"Training stopped early due to '{stopper_name}' condition.")


def log_training_complete(model_name: str, total_epochs: int) -> None:
    print(f"Training of '{model_name}' completed after {total_epochs} epochs.\n")


def print_training_summary(results: dict) -> None:
    print("-- -- Training Summary: -- --")
    print(f"Number of Epochs: {len(results['loss']['train'])}\n")

    print(f"Final Training Loss: {results['loss']['train'][-1]:.4f}")
    print(f"Final Validation Loss: {results['loss']['valid'][-1]:.4f}\n")

    print(f"Final Training Accuracy: {results['accuracy']['train'][-1]:.4f}")
    print(f"Final Validation Accuracy: {results['accuracy']['valid'][-1]:.4f}\n")


def is_training_mode(optimizer) -> bool:
    return bool(optimizer)


def get_context_manager(training_mode):
    return torch.set_grad_enabled(True) if training_mode else torch.inference_mode()


def calculate_average(metric, dataloader):
        return metric / len(dataloader)


def generate_evaluation_report(model, dataloader, best_training_step, parameters):
    evaluation_report = evaluate_model_performance(model, dataloader, "final")
    print_evaluation_report(evaluation_report)
    
    best_epoch = best_training_step["epoch"]

    if best_epoch < parameters["num_epochs"]:
        print(f"The best model weights from epoch {best_epoch} have been loaded.\n")
        model.load_state_dict(best_training_step["weights"])

        evaluation_report = evaluate_model_performance(model, dataloader, "best")
        print_evaluation_report(evaluation_report)
    
    return evaluation_report


def evaluate_model_performance(model, dataloader, evaluation_stage):
    loss_fn, _ = initialize_metrics()
    labels, predictions = fetch_labels_and_predictions(model, dataloader)

    report = prepare_classification_report(labels, predictions, loss_fn)
    report["evaluation_stage"] = evaluation_stage

    return report


def fetch_labels_and_predictions(model, dataloader):
    device = get_device_from_model(model)
    model.eval()

    all_labels, all_predictions = [], []

    with torch.inference_mode():
        for features, labels in dataloader:
            features, labels = transfer_data_to_device(features, labels, device)
            predictions = model(features)

            all_labels.append(labels)
            all_predictions.append(predictions)

    all_labels = torch.cat(all_labels, dim=0)
    all_predictions = torch.cat(all_predictions, dim=0)

    return all_labels, all_predictions


def prepare_classification_report(actual_labels, predictions, loss_fn):
    loss = loss_fn(predictions, actual_labels).item()
    predicted_labels = compute_predicted_labels(predictions)

    actual_labels = transfer_to_device(actual_labels, "cpu")
    predicted_labels = transfer_to_device(predicted_labels, "cpu")

    return {
        "loss": loss,
        "matrix": confusion_matrix(actual_labels, predicted_labels),
        "accuracy": accuracy_score(actual_labels, predicted_labels),
        "precision": precision_score(actual_labels, predicted_labels),
        "recall": recall_score(actual_labels, predicted_labels),
        "f1_score": f1_score(actual_labels, predicted_labels)
    }


def compute_predicted_labels(logits):
    label_probabilities = torch.softmax(logits, dim=1)
    predicted_labels = torch.argmax(label_probabilities, dim=1) 

    return predicted_labels


def print_evaluation_report(report):
    print("-- -- Model Evaluation -- --")
    print(f"-- Evaluation Stage: {report['evaluation_stage']} --")
    print(f"Loss: {report['loss']:.4f}")
    print(f"Accuracy: {report['accuracy']:.4f}\n")

    print(f"Confusion Matrix:\n {report['matrix']}\n")
    print(f"Precision Score: {report['precision']:.4f}")
    print(f"Recall Score: {report['recall']:.4f}")
    print(f"F1 Score: {report['f1_score']:.4f}\n")


def save_training_summary(target_dir, training_metrics, file_name="training_metrics.pkl"):
    with open(os.path.join(target_dir, file_name), "wb") as file:
        pickle.dump(training_metrics, file)


def save_evaluation_report(target_dir, evaluation_report, file_name="evaluation_report.pkl"):
    with open(os.path.join(target_dir, file_name), "wb") as file:
        pickle.dump(evaluation_report, file)
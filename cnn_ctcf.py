from setup import *

# Hypothetical Data
sample_sequences = ["ATCG", "CTAG", "CATG", "CGTA"]
sequence_labels = [1, 0, 1, 0]

# One-hot encoding
def convert_to_one_hot(sequence):
    conversion_map = {'A': [1, 0, 0, 0], 'C': [0, 1, 0, 0], 'G': [0, 0, 1, 0], 'T': [0, 0, 0, 1]}
    return [conversion_map[nucleotide.upper()] for nucleotide in sequence]
one_hot_sequences = [convert_to_one_hot(seq) for seq in sample_sequences]

# Splitting
X_train_set, X_test_set, y_train_set, y_test_set = train_test_split(one_hot_sequences, sequence_labels, test_size=0.25, random_state=42)
X_train_tensor = torch.tensor(X_train_set, dtype=torch.float32)
X_test_tensor = torch.tensor(X_test_set, dtype=torch.float32)
y_train_tensor = torch.tensor(y_train_set, dtype=torch.long)
y_test_tensor = torch.tensor(y_test_set, dtype=torch.long)


X_train_tensor = X_train_tensor.permute(0, 2, 1)
X_test_tensor = X_test_tensor.permute(0, 2, 1)

# CNN model
class SimpleCNN(nn.Module):
    def __init__(self):
        super(SimpleCNN, self).__init__()
        self.conv_layer1 = nn.Conv1d(in_channels=4, out_channels=64, kernel_size=1)
        self.fc_layer1 = nn.Linear(256, 2)

    def forward(self, x):
        x = self.conv_layer1(x)
        x = F.relu(x)
        x = x.view(x.size(0), -1)
        x = self.fc_layer1(x)
        return x

# Initialization, loss function, and optimizer
network_model = SimpleCNN()
loss_function = nn.CrossEntropyLoss()
optimizer_function = optim.Adam(network_model.parameters(), lr=0.001)

# Training
epoch_count = 30
for epoch_number in range(epoch_count):
    model_outputs = network_model(X_train_tensor)
    loss_value = loss_function(model_outputs, y_train_tensor)

    optimizer_function.zero_grad()
    loss_value.backward()
    optimizer_function.step()

    print(f'Epoch {epoch_number + 1}/{epoch_count}, Loss: {loss_value.item()}')

# Evaluation
network_model.eval()
with torch.no_grad():
    predictions = network_model(X_test_tensor)
    _, predicted_labels = torch.max(predictions, 1)
    model_accuracy = (predicted_labels == y_test_tensor).sum().item() / y_test_tensor.size(0)
    print(f'Model Accuracy: {model_accuracy:.2f}')

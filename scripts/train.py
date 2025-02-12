import torch
import torch.nn as nn
import torch.optim as optim
import sys
import os

# Add the project root directory to the Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


from app.models.model import SimpleModel  # Import your model class

# Create and train your model
def train_and_save_model():
    # Create sample dataset (replace with actual data if available)
    X = torch.randn(100, 10)  # 100 samples, 10 features each
    y = torch.randn(100, 1)   # 100 target values

    model = SimpleModel()
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=0.01)

    # Training loop
    for epoch in range(100):
        model.train()
        optimizer.zero_grad()
        output = model(X)
        loss = criterion(output, y)
        loss.backward()
        optimizer.step()

        if epoch % 10 == 0:
            print('Epoch {}, Loss: {}'.format(epoch, loss.item()))


    # Save the trained model's state
    torch.save(model.state_dict(), "app/models/trained_model.pth")
    print("Model saved at app/models/trained_model.pth")

if __name__ == "__main__":
    train_and_save_model()
